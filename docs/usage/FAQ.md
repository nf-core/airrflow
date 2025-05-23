# nf-core/airrflow: Frequently Asked Questions

## What are B cell or T cell clone?
B cell or T cell clone is a group of B cells or T cells that are derived from the same progenitor cell through clonal expansion. For T cells, this definition is more strict as T cells do not undergo somatic hypermutation, so the TCRs from T cells in the same clone should be identical. For B cells, on the other hand, the BCRs from cells in the same clone can differ due to somatic hypermutation. They also can have a variety of isotypes.


## What to consider for clonal analysis?

There are two crucial considerations for clonal analysis with nf-core/airrflow: : (1) the scope of sample grouping — that is, what samples should be grouped together before clonal inference; and (2) the clonal threshold, that is he maximum allowable distance between receptor sequences to be considered within the same clonal group. These considerations are explored in detail in the sections below.

### Sample grouping before defining clonal groups

We often want to analyze clonal groups from samples belonging to the same individual or animal model and being extracted across time, different conditions or from different tissues. To ensure that the same clone ID (field `clone_id` in the output AIRR rearrangement file) is assigned to the same BCR / TCR clone across these conditions, the clonal inference step should be done pulling the sequences from these samples together. This is why, by default, nf-core/airrflow uses the `subject_id` column to group samples prior to defining clonal groups, so it is important to set the exact same subject ID to samples from the same individual across different conditions.

The sample grouping can also be controlled with the [`--cloneby`](https://nf-co.re/airrflow/4.3.0/parameters/#cloneby) parameter, by providing the name of the column containing the group information that should be used to pull the samples together before defining clonal groups (samples or rows with the same string in this column will be grouped together). You can create a new column if you wish for this purpose.


### Setting a clonal threshold
nf-core/airrflow utilizes the Hierarchical clustering method in the Immcantation tool [SCOPer](https://scoper.readthedocs.io/) to infer clonal groups, which initially partitions the BCR / TCR sequences according to V gene, J gene and junction length. Then, it defines clonal groups within each partition by performing hierarchical clustering of the sequences within a partition and cutting the clusters according to an automatically detected or user-defined threshold. More details about this method can be found on the respective [SCOPer vignette](https://scoper.readthedocs.io/en/stable/vignettes/Scoper-Vignette/#).

By default, `--clonal_threshold` is set to be 'auto', allowing the clonal threshold to be determined automatically using a method included in the Immcantation tool [SHazaM](https://shazam.readthedocs.io/). You can read more details about the method in the [SHazaM vignette](https://shazam.readthedocs.io/en/stable/vignettes/DistToNearest-Vignette/).

The clonal threshold can also be customized through the `--clonal_threshold` parameter. It defines the maximum allowable distance between receptor sequences to be considered within the same clonal group. Specifically, this distance is measured as the length-normalized Hamming distance across the BCR junction regions.

For BCR data, we recommend using the default setting initially. After running the pipeline, you can review the automatically calculated threshold in the `find_threshold` report to make sure it fits the data appropriately. If the threshold is unsatisfactory, you can re-run the pipeline with a manually specified threshold (e.g. `--clonal_threshold 0.1`) that is appropriate for your data. For a low number of sequences that are insufficient to satisfactorily determine a threshold with this method, we generally recommend a threshold of 0.1 for human BCR data.

Since TCRs do not undergo somatic hypermutation, TCR clones are defined strictly by identical junction regions. For this reason, the `--clonal_threshold` parameter should be set to 0 for TCR data.


## How to include BCR lineage tree computation?

BCR lineage tree computation is performed using the Immcantation package [Dowser](https://dowser.readthedocs.io/). This step is skipped by default because it can be time-consuming depending on the size of the input data and the size of the clonal groups. To enable lineage tree computation, add the `--lineage_trees` parameter and set it to be `true`. You can easily add lineage tree computation to a previous analysis by re-running the pipeline with the `-resume` so all the previous analysis steps are cached and not recomputed.

Dowser supports different methods for the lineage tree computation, `raxml` is the default but you can set other methods with the `--lineage_tree_builder` parameter, and provide the software executable with the `--lineage_tree_exec` parameter.

## How to costumize the analysis and figures?

nf-core/airrflow is a standardized pipeline that performs the different computational analysis steps and provides standard figures for a first data exploration. The computations results (e.g. clonal inference, mutation frequency analysis) are stored in the output AIRR rearrangement repertoire files in newly generated columns under `clonal_analysis/define_clones/all_repertoires`. You can use these Airrflow results as input for customized analyses using R and the Immcantation tools. You can find the tutorial for Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html).


## How to update process resource requests?

By default the pipeline has set reasonable process resource requests (number of CPUs, RAM memory, time limits) to the compute system. Depending on the size of your datasets or your running infrastructure you can customize these requests. The `resourceLimits` option applies upper resource request limits to all the processes in the pipeline. Make sure to set resource limits that are not surpassing the available resources in your compute infrastructure. You can do so in the `resource.config` file provided with the `-c` parameter.

```json title="resource.config"
process {
   resourceLimits = [cpus: 8, memory: 72.GB, time: 24.h]
 }

To update the resource requests for a specific pipeline process, you can also provide specific process requests in this config file. For example, to update the resource requests for the `CHANGEO_ASSIGNGENES` process:

```json title="resource.config"
process {
   resourceLimits = [cpus: 8, memory: 72.GB, time: 24.h]

   withName:CHANGEO_ASSIGNGENES {
        cpus   = 2
        memory = 10.GB
        time   = 5h
   }
}
```

In nf-core pipelines, each process has a label indicating the resources that are being requested (`process_low`, `process_medium`, `process_high`, ...). The CPUs, RAM and time set up for each of these labels can be found in the [base.config](https://github.com/nf-core/airrflow/blob/master/conf/base.config) file. You can update the resource requests for all processes with a specific label by providing the updated configuration. For example here we update the process requests of processes with the `process_high` label:

```json title="resource.config"
process {
   resourceLimits = [cpus: 24, memory: 100.GB, time: 24.h]

   withLabel:process_high {
        cpus   = 24
        memory = 100.GB
        time   = 10h
   }
}
```

Note that the resource requests will never exceed what is specified in the `resourceLimits` line, so if you do want to increase the resource requests for specific processes, you should also increase the `resourceLimits` requests and run the pipeline in a compute infrastructure with sufficient resources. In this exmaple we also have updated the `resourceLimits` to reflect that.

> [!TIP]
> For more information about nf-core pipeline resource configurations, check out the [nf-core pipeline configuration docs](https://nf-co.re/docs/usage/getting_started/configuration).