import json
import sys
from pathlib import Path

import numpy as np
import tskit


def main() -> None:
    if len(sys.argv) != 2:
        raise SystemExit("Usage: python tskit-roundtrip.py OUTPUT_DIR")
    output_dir = Path(sys.argv[1])
    expected_count = int((output_dir / "case-count.txt").read_text().strip())
    tree_files = sorted(output_dir.glob("*.trees"))
    if len(tree_files) != expected_count:
        raise AssertionError((len(tree_files), expected_count))

    nested_mutations = 0
    for tree_file in tree_files:
        ts = tskit.load(tree_file)
        expected_file = tree_file.with_suffix(".tsv")
        if expected_file.exists():
            expected = np.loadtxt(
                expected_file, dtype=np.int8, delimiter="\t", ndmin=2
            )
            observed = ts.genotype_matrix().T
            if observed.shape != expected.shape or not np.array_equal(
                observed, expected
            ):
                mismatch = np.argwhere(observed != expected)
                raise AssertionError(
                    f"{tree_file.name}: observed={observed.shape}, "
                    f"expected={expected.shape}, mismatch={mismatch[:5]}"
                )
            if np.any(observed == tskit.MISSING_DATA):
                raise AssertionError(f"{tree_file.name}: missing sample genotype")
        elif ts.num_sites != 0 or ts.num_mutations != 0:
            raise AssertionError(f"{tree_file.name}: ancestry-only variants present")
        if ts.time_units != "generations":
            raise AssertionError(f"{tree_file.name}: unexpected time units")
        if [site.position for site in ts.sites()] != [
            locus + 0.5 for locus in range(ts.num_sites)
        ]:
            raise AssertionError(f"{tree_file.name}: unexpected site positions")

        expected_schema = {"codec": "json"}
        schemas = ts.table_metadata_schemas
        if (
            ts.metadata_schema.schema != expected_schema
            or schemas.individual.schema != expected_schema
            or schemas.node.schema != expected_schema
            or schemas.site.schema != expected_schema
        ):
            raise AssertionError(f"{tree_file.name}: JSON metadata schema missing")
        metadata = ts.metadata
        if metadata["software"]["name"] != "AlphaSimR":
            raise AssertionError(f"{tree_file.name}: invalid metadata")
        if metadata["coordinate_system"] != "locus_index":
            raise AssertionError(f"{tree_file.name}: invalid coordinate system")
        if (
            metadata["time_scale"] != "pedigree_depth"
            or metadata["time_origin"] != "deepest_recorded_individual"
            or metadata["founder_origin_offset_generations"] != 1
            or metadata["sample_ancestry_proxy_offset"]
            != "next_representable_older_time"
        ):
            raise AssertionError(f"{tree_file.name}: invalid time metadata")
        expected_encoding = (
            "founder_origins_and_sample_overrides"
            if expected_file.exists()
            else "none"
        )
        if metadata["variant_encoding"] != expected_encoding:
            raise AssertionError(f"{tree_file.name}: invalid variant encoding")
        for provenance in ts.provenances():
            record = json.loads(provenance.record)
            if record["schema_version"] != "1.0.0":
                raise AssertionError(f"{tree_file.name}: invalid provenance")
        for individual in ts.individuals():
            if not isinstance(individual.metadata, dict):
                raise AssertionError(f"{tree_file.name}: invalid individual metadata")
        for node in ts.nodes():
            node_metadata = node.metadata
            if node.is_sample() and node_metadata.get("role") != "sample":
                raise AssertionError(f"{tree_file.name}: invalid sample node")
        for locus, site in enumerate(ts.sites(), start=1):
            if site.metadata["locus"] != locus:
                raise AssertionError(f"{tree_file.name}: invalid site metadata")

        expected_samples = np.loadtxt(
            tree_file.with_suffix(".samples.tsv"),
            dtype=np.float64,
            delimiter="\t",
            ndmin=2,
        )
        observed_samples = []
        observed_ids = []
        for node_id in ts.samples():
            node = ts.node(node_id)
            node_metadata = node.metadata
            individual_metadata = ts.individual(node.individual).metadata
            observed_samples.append(
                [
                    individual_metadata["iid"],
                    node_metadata["homolog"],
                    node.time,
                ]
            )
            observed_ids.append(individual_metadata["id"])
        if not np.array_equal(np.asarray(observed_samples), expected_samples):
            raise AssertionError(
                f"{tree_file.name}: sample order or time changed"
            )
        expected_ids = [
            bytes.fromhex(line[1:]).decode("utf-8")
            for line in tree_file.with_suffix(".ids.txt").read_text().splitlines()
        ]
        if observed_ids != expected_ids:
            raise AssertionError(f"{tree_file.name}: sample IDs changed")

        is_full = (
            tree_file.stem.endswith("_full") or "_full_chr" in tree_file.stem
        )
        expected_pedigree = np.loadtxt(
            tree_file.with_suffix(".pedigree.tsv"),
            dtype=np.int64,
            delimiter="\t",
            ndmin=2,
        )
        if is_full:
            if ts.num_individuals != expected_pedigree.shape[0]:
                raise AssertionError(
                    f"{tree_file.name}: full pedigree was not retained"
                )
            observed_pedigree = []
            for expected_iid, individual in enumerate(ts.individuals(), start=1):
                individual_metadata = individual.metadata
                if individual_metadata["iid"] != expected_iid:
                    raise AssertionError(
                        f"{tree_file.name}: individual order changed"
                    )
                parents = [
                    0 if parent == tskit.NULL else parent + 1
                    for parent in individual.parents
                ]
                observed_pedigree.append(parents if parents else [0, 0])
            if not np.array_equal(
                np.asarray(observed_pedigree), expected_pedigree
            ):
                raise AssertionError(f"{tree_file.name}: pedigree changed")
        else:
            for individual in ts.individuals():
                individual_metadata = individual.metadata
                for parent in individual.parents:
                    if parent != tskit.NULL:
                        parent_metadata = ts.individual(parent).metadata
                        if parent_metadata["iid"] >= individual_metadata["iid"]:
                            raise AssertionError(
                                f"{tree_file.name}: invalid simplified pedigree"
                            )

        expected_origins = np.loadtxt(
            tree_file.with_suffix(".origins.tsv"),
            dtype=np.int64,
            delimiter="\t",
            ndmin=2,
        )
        observed_origins = np.empty_like(expected_origins)
        sample_nodes = ts.samples()
        for locus in range(expected_origins.shape[1]):
            tree = ts.at(locus + 0.5)
            for sample_index, sample_node in enumerate(sample_nodes):
                root = sample_node
                while tree.parent(root) != tskit.NULL:
                    root = tree.parent(root)
                root_metadata = ts.node(root).metadata
                observed_origins[sample_index, locus] = root_metadata[
                    "founder_haplotype_origin"
                ]
        if not np.array_equal(observed_origins, expected_origins):
            mismatch = np.argwhere(observed_origins != expected_origins)
            raise AssertionError(
                f"{tree_file.name}: ancestry mismatch={mismatch[:5]}"
            )
        nested_mutations += sum(
            mutation.parent != tskit.NULL for mutation in ts.mutations()
        )
        if tree_file.name.startswith(
            ("permuted_initial_founders_", "permuted_reset_founders_")
        ):
            sample_nodes = set(ts.samples())
            if any(mutation.node in sample_nodes for mutation in ts.mutations()):
                raise AssertionError(
                    f"{tree_file.name}: false terminal founder mutation"
                )

    if nested_mutations == 0:
        raise AssertionError("No recurrent/back mutation parent was validated")

    semantic_pairs = 0
    pairwise_mrcas = 0
    for simple_file in output_dir.glob("*_simple.trees"):
        full_file = Path(str(simple_file).replace("_simple.trees", "_full.trees"))
        if not full_file.exists():
            continue
        simple = tskit.load(simple_file)
        full = tskit.load(full_file)
        if (
            simple.num_samples != full.num_samples
            or simple.sequence_length != full.sequence_length
        ):
            raise AssertionError(
                f"{simple_file.name}: simplified dimensions changed"
            )
        simple_samples = simple.samples()
        full_samples = full.samples()
        simple_sample_index = {
            node: index for index, node in enumerate(simple_samples)
        }
        full_sample_index = {
            node: index for index, node in enumerate(full_samples)
        }
        breakpoints = sorted(
            set(simple.breakpoints()).union(full.breakpoints())
        )
        for left, right in zip(breakpoints, breakpoints[1:]):
            position = (left + right) / 2
            simple_tree = simple.at(position)
            full_tree = full.at(position)
            for first in range(simple.num_samples):
                for second in range(first + 1, simple.num_samples):
                    simple_mrca = simple_tree.mrca(
                        simple_samples[first], simple_samples[second]
                    )
                    full_mrca = full_tree.mrca(
                        full_samples[first], full_samples[second]
                    )
                    if simple_mrca == tskit.NULL:
                        simple_signature = None
                    else:
                        simple_signature = (
                            simple.node(simple_mrca).time,
                            tuple(
                                sorted(
                                    simple_sample_index[node]
                                    for node in simple_tree.samples(simple_mrca)
                                )
                            ),
                        )
                    if full_mrca == tskit.NULL:
                        full_signature = None
                    else:
                        full_signature = (
                            full.node(full_mrca).time,
                            tuple(
                                sorted(
                                    full_sample_index[node]
                                    for node in full_tree.samples(full_mrca)
                                )
                            ),
                        )
                    if simple_signature != full_signature:
                        raise AssertionError(
                            f"{simple_file.name}: simplification changed "
                            f"ancestry at {position} for samples "
                            f"{first}, {second}"
                        )
                    pairwise_mrcas += 1
        semantic_pairs += 1
    if semantic_pairs == 0 or pairwise_mrcas == 0:
        raise AssertionError("No simplification equivalence was validated")
    print(
        f"validated={len(tree_files)} tskit={tskit.__version__} "
        f"exact_alleles=True exact_ancestry=True sample_order=True "
        f"missing_genotypes=0 "
        f"nested_mutations={nested_mutations} "
        f"simplification_pairs={semantic_pairs} "
        f"pairwise_mrcas={pairwise_mrcas}"
    )


if __name__ == "__main__":
    main()
