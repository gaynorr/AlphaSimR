/**
Copyright Gary K. Chen (gchen98@gmail.com)

This file is part of the Markovian Coalescent Simulator (MaCS),
<https://github.com/gchen98/macs>. It has been modified by the AlphaSimR
authors for use in AlphaSimR. The file LICENSE.note, in the root of the
AlphaSimR sources, describes those modifications.

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

    http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
**/
#include<iostream>
using namespace std;

const string COMMAND="COMMAND:";
const string SEED="SEED:";
const string NEWICKTREE="NEWICK_TREE:";
const string HAPLOBEGIN="BEGIN_SEGREGATING_SITES";
const string HAPLOEND="END_SEGREGATING_SITES";
const string MUTATIONSITE="SITE:";
const string SNPBEGIN="BEGIN_SELECTED_SITES";
const string SNPEND="END_SELECTED_SITES";
const string TOTALSAMPLES="TOTAL_SAMPLES:";
const string TOTALSITES="TOTAL_SITES:";

const string FIELD_DELIMITER="\t";
const string TEMPFILE_PREFIX="haplo.";
const string TREEFILE_PREFIX="tree.";
