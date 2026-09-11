#include <RcppTskit.hpp>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <memory>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace {

void freeAlphaSimRTables(tsk_table_collection_t* tables) {
  if (tables != nullptr) {
    tsk_table_collection_free(tables);
    delete tables;
  }
}

using AlphaSimRTableXPtr =
    Rcpp::XPtr<tsk_table_collection_t, Rcpp::PreserveStorage,
               freeAlphaSimRTables, true>;

struct TableDeleter {
  void operator()(tsk_table_collection_t* tables) const {
    if (tables != nullptr) {
      tsk_table_collection_free(tables);
      delete tables;
    }
  }
};

void checkTsk(int ret) {
  if (ret < 0) {
    Rcpp::stop(tsk_strerror(ret));
  }
}

tsk_id_t checkTskId(tsk_id_t id) {
  if (id < 0) {
    Rcpp::stop(tsk_strerror(id));
  }
  return id;
}

std::string jsonEscape(const std::string& value) {
  std::ostringstream output;
  for (unsigned char character : value) {
    switch (character) {
      case '"': output << "\\\""; break;
      case '\\': output << "\\\\"; break;
      case '\b': output << "\\b"; break;
      case '\f': output << "\\f"; break;
      case '\n': output << "\\n"; break;
      case '\r': output << "\\r"; break;
      case '\t': output << "\\t"; break;
      default:
        if (character < 0x20) {
          output << "\\u" << std::hex << std::setw(4) << std::setfill('0')
                 << static_cast<int>(character) << std::dec;
        } else {
          output << character;
        }
    }
  }
  return output.str();
}

int getPloidy(SEXP history, int chromosome) {
  if (TYPEOF(history) == INTSXP) {
    return Rf_length(history);
  }
  if (TYPEOF(history) != VECSXP) {
    Rcpp::stop("Invalid recombination history entry");
  }
  Rcpp::List byChromosome(history);
  if (chromosome < 0 || chromosome >= byChromosome.size()) {
    Rcpp::stop("Chromosome is absent from recombination history");
  }
  SEXP chromosomeHistory = byChromosome[chromosome];
  if (TYPEOF(chromosomeHistory) != VECSXP) {
    Rcpp::stop("Invalid chromosome recombination history");
  }
  return Rf_length(chromosomeHistory);
}

}  // namespace

// [[Rcpp::export]]
Rcpp::List resolveVariantEncodingCpp(
    Rcpp::IntegerMatrix originRows,
    Rcpp::IntegerMatrix currentHaplotypes,
    Rcpp::IntegerMatrix originHaplotypes,
    Rcpp::LogicalVector knownOrigins) {
  const int nSamples = currentHaplotypes.nrow();
  const int nLoci = currentHaplotypes.ncol();
  const int nOrigins = originHaplotypes.nrow();
  if (originRows.nrow() != nSamples || originRows.ncol() != nLoci ||
      originHaplotypes.ncol() != nLoci ||
      knownOrigins.size() != nOrigins || nOrigins < 1) {
    Rcpp::stop("Variant-encoding dimensions disagree");
  }

  Rcpp::IntegerMatrix overrides(nSamples, nLoci);
  std::fill(overrides.begin(), overrides.end(), -1);
  std::vector<int> inferred(nOrigins, -1);
  for (int locus = 0; locus < nLoci; ++locus) {
    std::fill(inferred.begin(), inferred.end(), -1);
    for (int sample = 0; sample < nSamples; ++sample) {
      const int origin = originRows(sample, locus) - 1;
      const int allele = currentHaplotypes(sample, locus);
      if (origin < 0 || origin >= nOrigins) {
        Rcpp::stop("Sample ancestry references an unknown founder origin");
      }
      if (allele != 0 && allele != 1) {
        Rcpp::stop("Current haplotypes must contain only 0 and 1");
      }
      if (knownOrigins[origin] == NA_LOGICAL) {
        Rcpp::stop("Known-origin indicators cannot be missing");
      }
      if (!knownOrigins[origin]) {
        if (inferred[origin] == -1) {
          inferred[origin] = allele;
        } else if (inferred[origin] != allele) {
          inferred[origin] = -2;
        }
      }
    }
    for (int origin = 0; origin < nOrigins; ++origin) {
      if (!knownOrigins[origin] && inferred[origin] >= 0) {
        originHaplotypes(origin, locus) = inferred[origin];
      }
      const int allele = originHaplotypes(origin, locus);
      if (allele != 0 && allele != 1) {
        Rcpp::stop("Origin haplotypes must contain only 0 and 1");
      }
    }
    for (int sample = 0; sample < nSamples; ++sample) {
      const int origin = originRows(sample, locus) - 1;
      const int allele = currentHaplotypes(sample, locus);
      if (allele != originHaplotypes(origin, locus)) {
        overrides(sample, locus) = allele;
      }
    }
  }

  return Rcpp::List::create(
      Rcpp::Named("originHaplotypes") = originHaplotypes,
      Rcpp::Named("sampleAlleleOverrides") = overrides);
}

// [[Rcpp::export]]
SEXP buildTreeSequenceCpp(Rcpp::List recHist,
                          Rcpp::IntegerMatrix pedigree,
                          Rcpp::IntegerVector sampleIid,
                          Rcpp::CharacterVector individualId,
                          int chromosome,
                          int nLoci,
                          Rcpp::IntegerVector founderIid,
                          Rcpp::IntegerMatrix originHaplotypes,
                          Rcpp::IntegerMatrix sampleAlleleOverrides,
                          bool includeVariants,
                          bool simplify,
                          std::string version,
                          std::string timestamp) {
  const int nIndividuals = recHist.size();
  if (nIndividuals < 1 || pedigree.nrow() != nIndividuals ||
      pedigree.ncol() < 2 || individualId.size() != nIndividuals) {
    Rcpp::stop("Recombination history, pedigree, and individual IDs disagree");
  }
  if (nLoci < 1 || founderIid.size() < 1) {
    Rcpp::stop("Invalid tree-sequence dimensions");
  }

  std::vector<int> founders;
  std::unordered_set<int> founderSet;
  founders.reserve(founderIid.size());
  for (int iid : founderIid) {
    if (iid < 1 || iid > nIndividuals || !founderSet.insert(iid).second ||
        pedigree(iid - 1, 0) != 0 || pedigree(iid - 1, 1) != 0) {
      Rcpp::stop("Invalid founder identifiers");
    }
    founders.push_back(iid - 1);
  }

  std::unordered_set<int> samples;
  for (int iid : sampleIid) {
    if (iid < 1 || iid > nIndividuals) {
      Rcpp::stop("Sample iid is outside the recorded pedigree");
    }
    if (!samples.insert(iid - 1).second) {
      Rcpp::stop("Sample individual IDs must be unique");
    }
  }
  if (samples.empty()) {
    Rcpp::stop("At least one sample individual is required");
  }

  std::vector<int> depth(nIndividuals, 0);
  std::vector<int> ploidy(nIndividuals, 0);
  int maxDepth = 0;
  for (int i = 0; i < nIndividuals; ++i) {
    ploidy[i] = getPloidy(recHist[i], chromosome);
    if (ploidy[i] < 1) {
      Rcpp::stop("Recorded individual has zero ploidy");
    }
    const int mother = pedigree(i, 0);
    const int father = pedigree(i, 1);
    if (mother == 0 && father == 0) {
      depth[i] = 0;
    } else {
      if (mother < 1 || father < 1 || mother > i || father > i) {
        Rcpp::stop("Pedigree is not ordered before offspring");
      }
      depth[i] = std::max(depth[mother - 1], depth[father - 1]) + 1;
      maxDepth = std::max(maxDepth, depth[i]);
    }
  }

  std::vector<int> founderOrigins;
  std::unordered_map<int, int> originRow;
  for (int i : founders) {
    if (TYPEOF(recHist[i]) != INTSXP) {
      Rcpp::stop("Founder recombination history must contain haplotype IDs");
    }
    Rcpp::IntegerVector origins(recHist[i]);
    if (origins.size() != ploidy[i]) {
      Rcpp::stop("Founder haplotype IDs do not match founder ploidy");
    }
    for (int origin : origins) {
      if (origin < 1) {
        Rcpp::stop("Founder haplotype IDs must be positive");
      }
      if (originRow.emplace(origin, founderOrigins.size()).second) {
        founderOrigins.push_back(origin);
      }
    }
  }

  std::unique_ptr<tsk_table_collection_t, TableDeleter> tables(
      new tsk_table_collection_t());
  checkTsk(tsk_table_collection_init(tables.get(), 0));
  const std::string jsonSchema = "{\"codec\":\"json\"}";
  checkTsk(tsk_table_collection_set_metadata_schema(
      tables.get(), jsonSchema.c_str(), jsonSchema.size()));
  checkTsk(tsk_individual_table_set_metadata_schema(
      &tables->individuals, jsonSchema.c_str(), jsonSchema.size()));
  checkTsk(tsk_node_table_set_metadata_schema(
      &tables->nodes, jsonSchema.c_str(), jsonSchema.size()));
  checkTsk(tsk_site_table_set_metadata_schema(
      &tables->sites, jsonSchema.c_str(), jsonSchema.size()));
  tables->sequence_length = static_cast<double>(nLoci);
  const std::string timeUnits = "generations";
  checkTsk(tsk_table_collection_set_time_units(
      tables.get(), timeUnits.c_str(), timeUnits.size()));

  std::ostringstream metadata;
  metadata << "{\"software\":{\"name\":\"AlphaSimR\",\"version\":\""
           << jsonEscape(version) << "\"},\"chromosome\":" << chromosome + 1
           << ",\"coordinate_system\":\"locus_index\""
           << ",\"time_scale\":\"pedigree_depth\""
           << ",\"time_origin\":\"deepest_recorded_individual\""
           << ",\"founder_origin_offset_generations\":1"
           << ",\"sample_ancestry_proxy_offset\":"
           << "\"next_representable_older_time\""
           << ",\"variant_encoding\":\""
           << (includeVariants ? "founder_origins_and_sample_overrides" : "none")
           << "\"}";
  const std::string metadataString = metadata.str();
  checkTsk(tsk_table_collection_set_metadata(
      tables.get(), metadataString.c_str(), metadataString.size()));
  std::ostringstream provenance;
  provenance << "{\"schema_version\":\"1.0.0\",\"software\":{\"name\":"
             << "\"AlphaSimR\",\"version\":\"" << jsonEscape(version)
             << "\"},\"parameters\":{\"command\":\"asTreeSequence\""
             << ",\"chromosome\":" << chromosome + 1
             << ",\"coordinate_system\":\"locus_index\""
             << ",\"include_variants\":"
             << (includeVariants ? "true" : "false")
             << ",\"simplify\":" << (simplify ? "true" : "false") << "}}";
  const std::string provenanceString = provenance.str();
  checkTskId(tsk_provenance_table_add_row(
      &tables->provenances, timestamp.c_str(), timestamp.size(),
      provenanceString.c_str(), provenanceString.size()));

  std::vector<std::vector<tsk_id_t>> nodes(nIndividuals);
  std::vector<tsk_id_t> sampleNodes;
  std::vector<tsk_id_t> originNodes;
  originNodes.reserve(founderOrigins.size());
  for (std::size_t row = 0; row < founderOrigins.size(); ++row) {
    std::ostringstream originMetadata;
    originMetadata << "{\"founder_haplotype_origin\":"
                   << founderOrigins[row] << "}";
    const std::string originMetadataString = originMetadata.str();
    originNodes.push_back(checkTskId(tsk_node_table_add_row(
        &tables->nodes, 0, static_cast<double>(maxDepth + 1),
        TSK_NULL, TSK_NULL, originMetadataString.c_str(),
        originMetadataString.size())));
  }
  for (int i = 0; i < nIndividuals; ++i) {
    std::vector<tsk_id_t> parents;
    const int mother = pedigree(i, 0);
    const int father = pedigree(i, 1);
    if (mother > 0) {
      parents.push_back(mother - 1);
    }
    if (father > 0) {
      parents.push_back(father - 1);
    }
    const std::string label = Rcpp::as<std::string>(individualId[i]);
    std::ostringstream individualMetadata;
    individualMetadata << "{\"id\":\"" << jsonEscape(label)
                       << "\",\"iid\":" << i + 1 << "}";
    const std::string individualMetadataString = individualMetadata.str();
    const tsk_id_t individual = checkTskId(tsk_individual_table_add_row(
        &tables->individuals, 0, nullptr, 0,
        parents.empty() ? nullptr : parents.data(), parents.size(),
        individualMetadataString.c_str(), individualMetadataString.size()));
    if (individual != i) {
      Rcpp::stop("Unexpected tskit individual identifier");
    }
  }

  auto addIndividualNodes = [&](int i) {
    const bool isSample = samples.count(i) > 0;
    const double individualTime = static_cast<double>(maxDepth - depth[i]);
    const double lineageTime = isSample
        ? std::nextafter(individualTime,
                         std::numeric_limits<double>::infinity())
        : individualTime;
    nodes[i].reserve(ploidy[i]);
    for (int homolog = 0; homolog < ploidy[i]; ++homolog) {
      std::ostringstream nodeMetadata;
      nodeMetadata << "{\"homolog\":" << homolog + 1
                   << ",\"role\":\"ancestry\"}";
      const std::string nodeMetadataString = nodeMetadata.str();
      const tsk_id_t node = checkTskId(tsk_node_table_add_row(
          &tables->nodes, 0, lineageTime, TSK_NULL,
          isSample ? TSK_NULL : i, nodeMetadataString.c_str(),
          nodeMetadataString.size()));
      nodes[i].push_back(node);
    }
  };
  for (int i = 0; i < nIndividuals; ++i) {
    addIndividualNodes(i);
  }

  for (int iid : sampleIid) {
    const int i = iid - 1;
    const double sampleTime = static_cast<double>(maxDepth - depth[i]);
    for (int homolog = 0; homolog < ploidy[i]; ++homolog) {
      std::ostringstream nodeMetadata;
      nodeMetadata << "{\"homolog\":" << homolog + 1
                   << ",\"role\":\"sample\"}";
      const std::string nodeMetadataString = nodeMetadata.str();
      const tsk_id_t sampleNode = checkTskId(tsk_node_table_add_row(
          &tables->nodes, TSK_NODE_IS_SAMPLE, sampleTime, TSK_NULL, i,
          nodeMetadataString.c_str(), nodeMetadataString.size()));
      sampleNodes.push_back(sampleNode);
      checkTskId(tsk_edge_table_add_row(
          &tables->edges, 0, static_cast<double>(nLoci), nodes[i][homolog],
          sampleNode, nullptr, 0));
    }
  }

  for (int i : founders) {
    Rcpp::IntegerVector origins(recHist[i]);
    for (int homolog = 0; homolog < ploidy[i]; ++homolog) {
      const int row = originRow.at(origins[homolog]);
      checkTskId(tsk_edge_table_add_row(
          &tables->edges, 0, static_cast<double>(nLoci), originNodes[row],
          nodes[i][homolog], nullptr, 0));
    }
  }

  for (int i = 0; i < nIndividuals; ++i) {
    if (TYPEOF(recHist[i]) == INTSXP) {
      continue;
    }
    const int mother = pedigree(i, 0) - 1;
    const int father = pedigree(i, 1) - 1;
    if (mother < 0 || father < 0) {
      Rcpp::stop("Non-founder recombination record lacks two parents");
    }
    const int motherPloidy = ploidy[mother];
    const int fatherPloidy = ploidy[father];
    const int childPloidy = ploidy[i];
    int maternalCopies = 0;
    if (childPloidy == (motherPloidy + fatherPloidy) / 2) {
      maternalCopies = motherPloidy / 2;
    } else if (childPloidy == (motherPloidy + fatherPloidy) / 4) {
      maternalCopies = childPloidy;
    } else if (childPloidy == motherPloidy + fatherPloidy) {
      maternalCopies = motherPloidy;
    } else {
      Rcpp::stop("Unexpected parental ploidy levels");
    }

    Rcpp::List byChromosome(recHist[i]);
    Rcpp::List chromosomeHistory(byChromosome[chromosome]);
    for (int homolog = 0; homolog < childPloidy; ++homolog) {
      Rcpp::IntegerMatrix history(chromosomeHistory[homolog]);
      if (history.nrow() < 1 || history.ncol() != 2 || history(0, 1) != 1) {
        Rcpp::stop("Invalid recombination segment matrix");
      }
      const int parentIndividual = homolog < maternalCopies ? mother : father;
      for (int segment = 0; segment < history.nrow(); ++segment) {
        const int recordedHomolog = history(segment, 0);
        if (recordedHomolog >= 100) {
          Rcpp::stop(
              "Recombination history contains unresolved quadrivalent "
              "homolog labels and must be regenerated");
        }
        const int parentHomolog = recordedHomolog - 1;
        const int left = history(segment, 1) - 1;
        const int right = segment + 1 < history.nrow()
                              ? history(segment + 1, 1) - 1
                              : nLoci;
        if (parentHomolog < 0 ||
            parentHomolog >= static_cast<int>(nodes[parentIndividual].size()) ||
            left < 0 || right > nLoci || left >= right) {
          Rcpp::stop("Invalid parent homolog or recombination interval");
        }
        checkTskId(tsk_edge_table_add_row(
            &tables->edges, static_cast<double>(left),
            static_cast<double>(right), nodes[parentIndividual][parentHomolog],
            nodes[i][homolog], nullptr, 0));
      }
    }
  }

  if (includeVariants) {
    if (originHaplotypes.nrow() !=
            static_cast<int>(founderOrigins.size()) ||
        originHaplotypes.ncol() != nLoci) {
      Rcpp::stop("Origin haplotypes do not match recorded haplotype origins");
    }
    if (sampleAlleleOverrides.nrow() !=
            static_cast<int>(sampleNodes.size()) ||
        sampleAlleleOverrides.ncol() != nLoci) {
      Rcpp::stop("Sample allele overrides do not match sampled haplotypes");
    }
    for (int locus = 0; locus < nLoci; ++locus) {
      std::ostringstream siteMetadata;
      siteMetadata << "{\"locus\":" << locus + 1 << "}";
      const std::string siteMetadataString = siteMetadata.str();
      const tsk_id_t site = checkTskId(tsk_site_table_add_row(
          &tables->sites, static_cast<double>(locus) + 0.5, "0", 1,
          siteMetadataString.c_str(), siteMetadataString.size()));
      for (int row = 0; row < originHaplotypes.nrow(); ++row) {
        const int allele = originHaplotypes(row, locus);
        if (allele != 0 && allele != 1) {
          Rcpp::stop("Origin haplotypes must contain only 0 and 1");
        }
        if (allele == 1) {
          checkTskId(tsk_mutation_table_add_row(
              &tables->mutations, site, originNodes[row], TSK_NULL,
              TSK_UNKNOWN_TIME, "1", 1, nullptr, 0));
        }
      }
      for (int row = 0; row < sampleAlleleOverrides.nrow(); ++row) {
        const int allele = sampleAlleleOverrides(row, locus);
        if (allele < -1 || allele > 1) {
          Rcpp::stop("Sample allele overrides must contain only -1, 0, and 1");
        }
        if (allele >= 0) {
          const char* state = allele == 0 ? "0" : "1";
          checkTskId(tsk_mutation_table_add_row(
              &tables->mutations, site, sampleNodes[row], TSK_NULL,
              TSK_UNKNOWN_TIME, state, 1, nullptr, 0));
        }
      }
    }
  }

  checkTsk(tsk_table_collection_sort(tables.get(), nullptr, 0));
  checkTsk(tsk_table_collection_build_index(tables.get(), 0));
  if (includeVariants) {
    checkTsk(tsk_table_collection_compute_mutation_parents(tables.get(), 0));
  }
  checkTsk(tsk_table_collection_check_integrity(
      tables.get(), TSK_CHECK_TREES));
  if (simplify) {
    checkTsk(tsk_table_collection_simplify(
        tables.get(), sampleNodes.data(), sampleNodes.size(),
        TSK_SIMPLIFY_FILTER_INDIVIDUALS |
            TSK_SIMPLIFY_KEEP_INPUT_ROOTS,
        nullptr));
    checkTsk(tsk_table_collection_sort(tables.get(), nullptr, 0));
    checkTsk(tsk_table_collection_build_index(tables.get(), 0));
    checkTsk(tsk_table_collection_check_integrity(
        tables.get(), TSK_CHECK_TREES));
  }

  AlphaSimRTableXPtr xptr(tables.release(), true);
  return xptr;
}
