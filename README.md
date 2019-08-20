# rpRanker

## Docker Base image

Create the base image for all the RetroPath2.0 tools using the rpRanker scripts

```
docker build -t ibisba/rpbase -f Dockerfile .
```

## Requirements

* [libsbml](http://sbml.org/Software/libSBML) - SBML reader/writer
* [rdkit](https://www.rdkit.org/) - Chemoinformatics software
* [openbabel](http://openbabel.org/wiki/Main_Page) - Chemical toolbox
* [pycobra](https://modal.lille.inria.fr/pycobra/) - FBA simulation toolbox
* [scipy](https://www.scipy.org/) - Scientific package

## Native run

NOTE: one must import rpRanker on the root of the rpRanker folder

### Cache

The ranker requires pre-parsing a number of files (listed in input_cache README.md) from equilibrator, retrorules and MetaNetX by running:

```
python rpCache.py
```

It generates a folder with a collection of pickle files used by the program.

### Parse and Generate SBML's

To analyse the output of RetroPath2.0 (combined with RP2paths) and RetroPath3.0, first parse the outputs to SBML's using rpReader:

```
import rpRanker
rpreader = rpRanker.rpReader()
```

Then depending on your RP version follow the instructions:

#### RetroPath2.0

```
rpreader.compounds(rp2paths_compounds_file)
rpreader.transformation(rp2paths_scope)
rpreader.outPaths(rp2paths_out_paths, max_rule_ids)
rpreader.pathsToSBML()
rpreader.sbml_paths
```

max_rule_ids (default 10) represents the maximal allowed subpaths to be returned. Indeed RetroPath2.0 may have multiple rules associated with each reaction. rpReader takes that into consideration by enumerating all the possible pathways that may be formed from the possible different cofactors associated with the reactions rules. To avoid a combinatorial explosion, the rules with the highest score are used. The results are contained within the rpreader.sbml_paths parameter.

#### RetroPath3.0 

For a single JSON file:

```
rpsbml = rpreader.jsonToSBML(json_dict)
```

The input must be a parsed JSON as a dictionnary (not a file). It returns a rpSBML object.

For RetroPath3.0 with a collection of JSON's (batch mode):

```
rpreader.collectionJSON(collection_json_dict)
```

The input is a dictionnary of JSON dictionnaries and the results are contained within the rpreader.sbml_paths parameter.

To export the results as a TAR one can use the tools package:

```
tools = rpRanker.tools()
tools.writerpSBMLtar(rpreader.sbml_paths, path_to_tar)
``` 

### Add the cofactors

Open the collections of rpSBML's as a dictionnary:

```
rptools = rpRanker.tools()
rpsbml_paths = rptools.readrpSBMLtar(path_to_tar)
``` 

Or pass your own dictionnary of rpSBML objects and add the cofactors:

```
rpcofactors = rpRanker.rpCofactors()
for rpsbml_name in rpsbml_paths:
    rpcofactors.addCofactors(rpsbml.paths[rpsbml_name])
```

To export the results as a TAR one can use the tools package:

```
rptools.writerpSBMLtar(rpsbml_paths, path_to_tar)
``` 

### Thermodynamics

Open the collections of rpSBML's as a dictionnary:

```
rptools = rpRanker.tools()
rpsbml_paths = rptools.readrpSBMLtar(path_to_tar)
``` 

Or pass your own dictionnary of rpSBML objects and calculate the thermodynamics:

```
rpthermo = rpRanker.rpThermo()
for rpsbml_name in rpsbml_paths:
    rpthermo.pathway_drG_prime(rpsbml.paths[rpsbml_name])
```

To export the results as a TAR one can use the tools package:

```
rptools.writerpSBMLtar(rpsbml_paths, path_to_tar)
``` 

### Selenzyme

Open the collections of rpSBML's as a dictionnary:

```
rptools = rpRanker.tools()
rpsbml_paths = readrpSBMLtar(path_to_tar)
``` 

Or pass your own dictionnary of rpSBML objects, then query the Selenzyme web tool:

```
for rpsbml_name in rpsbml_paths:
    rptools.rpSelenzyme(rpreader.sbml_paths[rpsbml_name], selenzyme_rest_url)
```

To export the results as a TAR one can use the tools package:

```
tools.writerpSBMLtar(rpsbml_paths, path_to_tar)
``` 

### FBA

Open the collections of rpSBML's as a dictionnary:

```
rptools = rpRanker.tools()
rpsbml_paths = readrpSBMLtar(path_to_tar)
``` 

Or pass your own dictionnary of rpSBML objects, then calculate the flux using:

```
for rpsbml_name in rpsbml_paths:
    rptools.runFBA(rpsbml_paths[rpsbml_name], path_to_gem_model) 
```

gen_model is a genome scale model (GEM) but one can pass any SBML file to it. However, for it to work one needs to have an SBML model that include the same chemical species as the sink provided to RP2 or RP3.

To export the results as a TAR one can use the tools package:

```
rptools.writerpSBMLtar(rpsbml_paths, path_to_tar)
``` 

### Report

Open the collections of rpSBML's as a dictionnary:

```
tools.writeReport(path_to_tar, path_to_csv)
```

### Generate Sink

Generate the sink (for RetroPath2.0) from a GEM model:

```
rpsbml = rpSBML('tmp')
rpsbml.readSBML(params.inSBML)
rpreader = rpRanker.rpReader()
rpreader.genSink(rpsbml, path_to_csv, compartmentId)
```

path_to_csv is the path to the output CSV sink file, while compartmentId is the compartment from which the sink will be generated.

### Mathilde

```
import rpRanker
json_data = open(path_to_json)
json_dict = json.load(json_data)
rpreader = rpRanker.rpReader()
rpsbml = rpreader.jsonToSBML(json_dict)
rpcofactors = rpRanker.rpCofactors()
rpcofactors.addCofactors(rpsbml)
rptools = rpRanker.tools()
rptools.runFBA(rpsbml, path_to_gem_model)
target_reaction = rpsbml.model.getReaction('targetSink')
ibisba_annot_dict = rpsbml.readIBISBAAnnotation(target_reaction.getAnnotation())
ibisba_annot_dict['fba_rpFBA_obj'] <-- thats the RP target flux
```

fba_rpFBA_obj is the RP target flux using it as an objective. However, if the GEM model contains an objective using the FBC package or indeed if another objective was manually added, one can access the fluxes of the target under different objective function by selecting its name.

