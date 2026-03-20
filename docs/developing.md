# Development documentation

## Pipeline specific conventions

- [ ] Use `nextflow` autoformatter
- [ ] Channel assingment with `=` NOT `.set`
- [ ] Parameter names structure: `<pipelinesection>_<toolname>_<parameter>`
  - All components of structure should not have `_`, i.e., they should all be concatenated: `amp_hmmsearch_savealignments` not `amp_hmmsearch_save_alignments`
  - Exception: subworkflow specific activation parameters must start with 'run': `run_<pipelinesection>_screening`, e.g. `run_arg_screening`
  - Execption: tool specific skipping parameters must use 'skip' in the second position: `<pipelinesection>_skip_<toolname>`, e.g. `amp_skip_macrel`

## Adding new tool workflow

This checklist covers adding a specific tool to an _existing_ screening subworkflow.

> [!NOTE]
> Does not have to be in this precise order

- [ ] Installed modules `nf-core modules install <tool>/<subtool>`
- [ ] Added tools(s) to relevant `subworkflows/local/<screentype>.nf`
  - [ ] Added relevant modules at the top using `include` statement
  - [ ] Added the tool-specific if/else statement controlled with `params.<screentype>_skip_<toolname>`
  - [ ] (Old nf-core module template only) Version channels mixed for every newly module
  - [ ] (If applicable) Include auto-database downloading module, if tool needs it
  - [ ] (If applicable) Within if/else condition, output channels mixed into the subworkflow's final aggregation/summary tool mix output file for aggregating tool (e.g. hAMRronize, AMPcombi etc.)
  - [ ] (If applicable) MultiQC channels mixed
- [ ] Updated `workflows/funcscan.nf`
  - [ ] (If applicable) Added new input files via dedicated input channel
  - [ ] (If applicable) Added tool specific input control conditions within the screening subworkflow's if/else statement
- [ ] Added parameters and defaults added to `nextflow.config`
  - [ ] Added all new parameters following [pipeline-specific conventions](pipeline-specific conventions)
  - [ ] (If applicable) Where possible, include parameter for supplying locally-downloaded database
- [ ] Update `modules.conf`
  - [ ] Added `withName:` block
  - [ ] Added `ext.args = {}` entry including relevant new pipeline parameters
  - [ ] Added `publishDir` placing relevant output files into the screening-subworkflow specific directory with `"${params.outdir}/<screentype>/<toolname>/`
- [ ] If necessary, added any additional pre-execution parameter validation checks at the top of `subworkflows/local/utils_nfcore_funcscan_pipeline.nf` (e.g. for mutually exclusive parameters)
- [ ] Updated Documentation
  - [ ] `nf-core pipelines schema build` has been run and updated
    - [ ] Checked all tool-specific pipeline parameters moved to the relevant screen type section of schema
    - [ ] Checked all tool-specific pipeline parameters have short help text
    - [ ] (If appplicable) Checked all tool-specific pipeline parameters have validation checks added (e.g. number range, fixed list etc.)
    - [ ] Checked all tool-specific pipeline parameters have long-description with more information, including pointing to original documentation of tool itself
    - [ ] (If applicable) Checked all tool-specific pipeline parameters have the `Modifies tool parameter(s)` quote block
  - [ ] Added citation to `CITATIONS.md` (citation style: APA 7th edition)
  - [ ] Added citation to the toolCitation/BibliographyText functions in `subworkflows/local/utils_nfcore_funcscan_pipeline`
    - [ ] Added in-text citation
    - [ ] Added bibliography (citation style: APA 7th edition)
  - [ ] Added relevant documentation to `usage.md`
    - [ ] (If applicable) Added entry in 'Databases and reference files' on how to download databases manually
    - [ ] (If applicable) Added entry in 'Notes on screening tools <...>' if specific guidance is needed for execution
    - [ ] (If applicable) If new input sample input files (e.g. annotation files) required, updated samplesheet description
  - [ ] Described module output in `output.md`
    - [ ] Added entry in relevant screening subworkflow section in introduction
    - [ ] Added entry in introduction `tree` of whole output directory
    - [ ] Added entry in 'Pipeline overview' table of contents
    - [ ] Added dedicated 'Tool details' section in relevant subworkflow section including collapsable output list, description of a tool, and (ideally) description what primary output files can be used for
      - [ ] Checked all output files specified in the `pattern:` section of `publishDir` are listed if `pattern:` is used, otherwise just all files found in results directory
  - [ ] Added entry to 'Pipeline summary list' on `README.md`
  - [ ] (Optional) Added to pipeline metro map diagram (can be done just before release)
  - [ ] (First time contributor) add or move yourself to the Team list on `README.md`!
  - [ ] (First time contributor) add or move yourself to the manifest section of `nextflow.config` as `contributor`
- [ ] (If applicable) On nf-core/test-data: added small test-database on the [funcscan](https://github.com/nf-core/test-datasets/tree/funcscan) branch
  - [ ] Added documentation of source and/or how test-data generated/modified
- Updated relevant `conf/test*.conf`
  - Specified skipping new tool to `true` in `test_minimal.config`
  - (If applicable) Included/adjusted parameters in all relevant test configs
  - (If applicable) Added paths to new test database files
- Updated tests
  - Added assertions to all output files for each test the tool is executed in
  - Updated relevant snapshots with `nf-test --tag <tag> --profile +<docker,conda/apptainer> --update-snapshot`
  - Checked assertions stable with `nf-test --tag <tag> --profile +<docker,conda/apptainer>`
- Added entry to `CHANGELOG.md` (note: PR number can be added after)
  - Tagged issue reporter/feature requester as well as author of PR

## Adding a new screening subworkflow workflow

Screening subworkflows should

- [ ] Be written as a local subworkflow
- [ ] By default run all tools, and offer skipping of execution
- [ ] Include final emit channels at a minimum including:
  - [ ] Versions channel
  - [ ] (If any tools supported) MultiQC channel
  - [ ] (If aggregation tool exists) include an aggregation tool step

Subworkflows within the primary `workflow/funcscan.nf` file, should

- [ ] Should be imported at the top of the module
- [ ] Have a dedicated if/else statement with running with and without taxonomic classification
  - [ ] Should include a empty file filter
  - [ ] (If subworkflow includes tools using old nf-core/modules structure) Should include the versions mixing
- [ ] (If applicable) be added to the Annotation if/else statement
- [ ] (If applicable) be included in the MultiQC annotation if/else statement
