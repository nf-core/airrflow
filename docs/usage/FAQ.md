# nf-core/airrflow: Frequently Asked Questions

## How to update process resource requests and resource limits?

By default, the pipeline defines reasonable resource requests for each process (number of CPUs, RAM memory, time limits) based on typical compute environments. However, you can adjust these settings to better match the size of your datasets or the capabilities of your compute infrastructure. You can customize the limits and requests in `resource.config` file and provide it to the pipeline using the -c parameter during execution. The `resourceLimits` option applies upper resource request limits to all the processes in the pipeline. Ensure that these limits do not exceed the available resources on your compute system.

```json title="resource.config"
process {
   resourceLimits = [cpus: 8, memory: 72.GB, time: 24.h]
 }
```

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

Note that the resource requests will never exceed what is specified in the `resourceLimits` line, so if you do want to increase the resource requests for specific processes, you should also increase the `resourceLimits` requests and run the pipeline in a compute infrastructure with sufficient resources. In this example we also have updated the `resourceLimits` to reflect that.

> [!TIP]
> For more information about nf-core pipeline resource configurations, check out the [nf-core pipeline configuration docs](https://nf-co.re/docs/usage/getting_started/configuration).

## How to customize the analysis and figures?

nf-core/airrflow is a standardized pipeline that performs the different computational analysis steps and provides standard figures for a first data exploration. You can use these Airrflow results as input for customized analyses using R and the Immcantation tools. You can find the tutorial for Immcantation's single-cell V(D)J analysis [here](https://immcantation.readthedocs.io/en/stable/getting_started/10x_tutorial.html).
