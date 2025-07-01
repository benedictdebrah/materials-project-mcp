Materials Project API 0.84.3rc4.dev9+g7e0ad853 

OAS 3.1
--------------------------------------------------------

[/openapi.json](/openapi.json)

The Materials Project API allows anyone to have direct access to current, up-to-date information from the Materials Project database in a structured way.

This allows for analysis, development of automated tools, machine learning, downloading personal copies of the Materials Project database and more on a large scale.

The API is offered with the hopes of making Materials Project data more useful to you. We want you to use our data! As such, the API is offered free-of-charge and we support several tools to help you get started.

API Key
-------

To make any request to the Materials Project API, you must use an API key. Your API key is generated for you automatically upon registering with the Materials Project website and is synced with the email you used to register.

Remember to keep your API key safe and to not share it with anyone you do not trust.

If you are logged in, you can always access your API key from this page or from your [dashboard](https://next-gen.materialsproject.org/dashboard)
.

If you intend heavy API usage, you can give us a heads up by sending a message to [heavy.api.use@materialsproject.org](mailto:heavy.api.use@materialsproject.org)
. With the exception of retrieving charge densities, this is not required, but may help us if we see unusual load on our servers.

Accessing Data
--------------

To use the API, you have three options:

1.  You can use our first-party supported Python client. This is the recommend route. The `mp-api` package containing the client is pip installable.
    
        pip install mp-api
        
    
    The `MPRester` client can be accessed by importing from it. This will ultimately replace the legacy `MPRester` available in pymatgen.
    
    For more details on how to use this, including code examples, please see [https://next-gen.materialsproject.org/api](https://next-gen.materialsproject.org/api)
    .
    
2.  You can demo the API interactively on this documentation page. Click the "Authorize" button, paste in your API key, and then click the appropriate section to try out a query.
    
3.  Since this is a REST API, and offers a fully-compliant OpenAPI specification, it's possible to use the API with many libraries in many languages and environments, including JavaScript, MATLAB, Mathematica, etc. However, we do not offer first-party support for explaining how to do this, and you will have to follow the specification yourself.
    

Authorize

### [Materials Summary](#/Materials%20Summary)

Route providing a large amount of amalgamated data for a material. This is constructed by combining subsets of data from many of the other API endpoints. The summary endpoint is very useful for performing queries for materials over a large property space. Note that every unique material within the Materials Project should have a set of summary data. See the `SummaryDoc` schema for a full list of fields returned by this route.

GET

[/materials/summary/stats/](#/Materials%20Summary/search_materials_summary_stats__get)

Get SummaryStats documents

GET

[/materials/summary/](#/Materials%20Summary/search_materials_summary__get)

Get SummaryDoc documents

### [Materials](#/Materials)

Route for "core" information associated with a given material in the Materials Project database. The unique identifier for a material is its `material_id` (e.g. `mp-149`). Core data in this context refers to the crystal structure, information associated with it such as the density and chemical formula, and the associated calculations which are identified with unique `task_id` values. It does not contain any materials properties such as the formation energy or band gap, please consult other property-specific endpoints for this information. See the `MaterialsDoc` schema for a full list of fields returned by this route.

GET

[/materials/core/blessed\_tasks/](#/Materials/search_materials_core_blessed_tasks__get)

Get MaterialsDoc documents

POST

[/materials/core/find\_structure/](#/Materials/search_materials_core_find_structure__post)

Post FindStructure documents

GET

[/materials/core/formula\_autocomplete/](#/Materials/search_materials_core_formula_autocomplete__get)

Get FormulaAutocomplete documents

GET

[/materials/core/](#/Materials/search_materials_core__get)

Get MaterialsDoc documents

### [Materials Tasks](#/Materials%20Tasks)

Route for "core" information associated with a given calculation in the Materials Project database. Multiple calculations can ultimately be associated with a unique material, and are the source of its reported properties. The unique identifier for a calculation is its `task_id`. Note that the `material_id` chosen for a given material is sourced from one of the `task_id` values associated with it. Core data in this context refers to calculation quantities such as parsed input and output data (e.g. VASP input flags, atomic forces, structures) and runtime statistics. See the `TaskDoc` schema for a full list of fields returned by this route.

GET

[/materials/tasks/trajectory/](#/Materials%20Tasks/search_materials_tasks_trajectory__get)

Get TrajectoryDoc documents

GET

[/materials/tasks/entries/](#/Materials%20Tasks/search_materials_tasks_entries__get)

Get EntryDoc documents

GET

[/materials/tasks/deprecation/](#/Materials%20Tasks/search_materials_tasks_deprecation__get)

Get DeprecationDoc documents

GET

[/materials/tasks/](#/Materials%20Tasks/search_materials_tasks__get)

Get TaskDoc documents

### [Materials Thermo](#/Materials%20Thermo)

Route providing computed thermodynamic data for a material such as formation energy and energy above hull. Corrected energy values are also available that employ the schemes discussed by [Jain _et al._](https://journals.aps.org/prb/abstract/10.1103/PhysRevB.84.045115)
 and [Wang _et al._](https://chemrxiv.org/engage/chemrxiv/article-details/60c758d9469df42a4ef45757)
 See the `ThermoDoc` schema for a full list of fields returned by this route.

[For a more detailed description](https://docs.materialsproject.org/methodology/total-energies)

GET

[/materials/thermo/phase\_diagram/](#/Materials%20Thermo/search_materials_thermo_phase_diagram__get)

Get PhaseDiagramDoc documents

GET

[/materials/thermo/](#/Materials%20Thermo/search_materials_thermo__get)

Get ThermoDoc documents

### [Materials Dielectric](#/Materials%20Dielectric)

Route providing computed dielectric data for a material following the methodology discussed by [Petousis _et al._](https://doi.org/10.1038/sdata.2016.134)
 Note that dielectric data has not been calculated for all materials in the Materials Project database. See the `DielectricDoc` schema for a full list of fields returned by this route.

[For a more detailed description](https://docs.materialsproject.org/methodology/dielectricity)

GET

[/materials/dielectric/](#/Materials%20Dielectric/search_materials_dielectric__get)

Get DielectricDoc documents

### [Materials Magnetism](#/Materials%20Magnetism)

Route providing computed magnetic ordering related data for a material following the methodology discussed by [Horton _et al._](https://doi.org/10.1038/s41524-019-0199-7)
 Note that magnetic data has not been calculated for all materials in the Materials Project database. See the `MagnetismDoc` schema for a full list of fields returned by this route.

GET

[/materials/magnetism/](#/Materials%20Magnetism/search_materials_magnetism__get)

Get MagnetismDoc documents

### [Materials Piezoelectric](#/Materials%20Piezoelectric)

Route providing computed piezoelectric data for a material following the methodology discussed by [de Jong _et al._](https://doi.org/10.1038/sdata.2015.53)
 Note that piezoelectric data has not been calculated for all materials in the Materials Project database. See the `PiezoDoc` schema for a full list of fields returned by this route.

[For a more detailed description](https://docs.materialsproject.org/methodology/piezoelectricity)

GET

[/materials/piezoelectric/](#/Materials%20Piezoelectric/search_materials_piezoelectric__get)

Get PiezoelectricDoc documents

### [Materials Phonon](#/Materials%20Phonon)

Route providing computed phonon data for a material following the methodology discussed by [Petretto _et al._](https://doi.org/10.1038/sdata.2018.65)
 Note that phonon data has not been calculated for all materials in the Materials Project database. See the `PhononBSDOSDoc` schema for a full list of fields returned by this route.

[For a more detailed description](https://docs.materialsproject.org/methodology/phonons)

GET

[/materials/phonon/](#/Materials%20Phonon/search_materials_phonon__get)

Get PhononBSDOSDoc documents

### [Materials EOS](#/Materials%20EOS)

Route providing computed equations of state data for a material following the methodology discussed by [Latimer _et al._](https://doi.org/10.1038/s41524-018-0091-x)
 Note that equations of state data has not been calculated for all materials in the Materials Project database. See the `EOSDoc` schema for a full list of fields returned by this route.

[For a more detailed description](https://docs.materialsproject.org/methodology/equations-of-state)

GET

[/materials/eos/](#/Materials%20EOS/search_materials_eos__get)

Get EOSDoc documents

### [Materials Similarity](#/Materials%20Similarity)

Route providing a computed similarity metric between materials following the methodology discussed by Zimmerman _et al._ in [10.3389/fmats.2017.00034](https://doi.org/10.3389/fmats.2017.00034)
 and [10.1039/C9RA07755C](https://doi.org/10.1039/C9RA07755C)
. Note that similarity data has not been calculated for all materials in the Materials Project database. See the `imilarityDoc` schema for a full list of fields returned by this route.

GET

[/materials/similarity/](#/Materials%20Similarity/search_materials_similarity__get)

Get SimilarityDoc documents

### [Materials XAS](#/Materials%20XAS)

Route providing computed x-ray absorption spectroscopy data for a material following the methodology discussed by [Mathew _et al._](https://doi.org/10.1038/sdata.2018.151)
 and [Chen _et al._](https://doi.org/10.1038/s41597-021-00936-5)
 Note that x-ray absorption spectroscopy data has not been calculated for all materials in the Materials Project database. See the `XASDoc` schema for a full list of fields returned by this route.

GET

[/materials/xas/](#/Materials%20XAS/search_materials_xas__get)

Get XASDoc documents

### [Materials Grain Boundaries](#/Materials%20Grain%20Boundaries)

Route providing computed grain boundary data for a material following the methodology discussed by [Hui _et al._](https://doi.org/10.1016/j.actamat.2019.12.030)
 Note that grain boundary data has not been calculated for all materials in the Materials Project database. See the `GrainBoundaryDoc` schema for a full list of fields returned by this route.

GET

[/materials/grain\_boundaries/](#/Materials%20Grain%20Boundaries/search_materials_grain_boundaries__get)

Get GrainBoundaryDoc documents

### [Materials Electronic Structure](#/Materials%20Electronic%20Structure)

Routes providing computed electronic structure related data for a material such as band gap and fermi level. Python objects for line-mode band structures, density of states, and fermi surfaces are also available. This data was obtained following the methodology discussed by [Munro _et al._](https://doi.org/10.1038/s41524-020-00383-7)
 and [Ganose _et al._](https://doi.org/10.21105/joss.03089)
 Note that full band structure, density of states, and fermi surface data has not been calculated for all materials in the Materials Project database. See the `ElectronicStructureDoc` and `FermiDoc` schema for a full list of fields returned by the associated routes.

[For a more detailed description](https://docs.materialsproject.org/methodology/electronic-structure)

GET

[/materials/fermi/](#/Materials%20Electronic%20Structure/search_materials_fermi__get)

Get FermiDoc documents

GET

[/materials/electronic\_structure/bandstructure/](#/Materials%20Electronic%20Structure/search_materials_electronic_structure_bandstructure__get)

Get ElectronicStructureDoc documents

GET

[/materials/electronic\_structure/dos/](#/Materials%20Electronic%20Structure/search_materials_electronic_structure_dos__get)

Get ElectronicStructureDoc documents

GET

[/materials/electronic\_structure/](#/Materials%20Electronic%20Structure/search_materials_electronic_structure__get)

Get ElectronicStructureDoc documents

GET

[/materials/electronic\_structure/bandstructure/object/](#/Materials%20Electronic%20Structure/search_materials_electronic_structure_bandstructure_object__get)

Get BSObjectDoc documents

GET

[/materials/electronic\_structure/dos/object/](#/Materials%20Electronic%20Structure/search_materials_electronic_structure_dos_object__get)

Get DOSObjectDoc documents

### [Materials Elasticity](#/Materials%20Elasticity)

Route providing computed elasticity data for a material following the methodology discussed by [de Jong _et al._](https://doi.org/10.1038/sdata.2015.9)
 Note that elasticity data has not been calculated for all materials in the Materials Project database. See the `ElasticityDoc` schema for a full list of fields returned by this route.

[For a more detailed description](https://docs.materialsproject.org/methodology/elasticity)

GET

[/materials/elasticity/](#/Materials%20Elasticity/search_materials_elasticity__get)

Get ElasticityDoc documents

### [Materials Substrates](#/Materials%20Substrates)

Route providing computed suggested substrate data for a material following the methodology discussed by [Ding _et al._](https://doi.org/10.1021/acsami.6b01630)
 Note that substrate data has not been calculated for all materials in the Materials Project database. See the `SubstratesDoc` schema for a full list of fields returned by this route.

GET

[/materials/substrates/](#/Materials%20Substrates/search_materials_substrates__get)

Get SubstratesDoc documents

### [Materials Surface Properties](#/Materials%20Surface%20Properties)

Route providing computed surface property data for a material following the methodology discussed by [Tran _et al._](https://doi.org/10.1038/sdata.2016.80)
 Note that surface data has not been calculated for all materials in the Materials Project database. See the `SurfacePropDoc` schema for a full list of fields returned by this route.

GET

[/materials/surface\_properties/](#/Materials%20Surface%20Properties/search_materials_surface_properties__get)

Get SurfacePropDoc documents

### [Materials Robocrystallographer](#/Materials%20Robocrystallographer)

Route providing a computed text description for a material following the methodology discussed by [Ganose _et al._](https://doi.org/10.1557/mrc.2019.94)
 Note that descriptions may not been calculated for all materials in the Materials Project database. See the `RobocrysDoc` schema for a full list of fields returned by this route.

GET

[/materials/robocrys/text\_search/](#/Materials%20Robocrystallographer/search_materials_robocrys_text_search__get)

Get RobocrystallogapherDoc documents

GET

[/materials/robocrys/](#/Materials%20Robocrystallographer/search_materials_robocrys__get)

Get RobocrystallogapherDoc documents

### [Materials Synthesis](#/Materials%20Synthesis)

Route providing a synthesis recipes for materials extracted from literature following the methodology discussed by [Kononova _et al._](https://doi.org/10.1038/s41597-019-0224-1)
 Note that synthesis recipes may not be available for all materials in the Materials Project database. See the `SynthesisSearchResultModel` schema for a full list of fields returned by this route.

GET

[/materials/synthesis/](#/Materials%20Synthesis/search_materials_synthesis__get)

Get SynthesisSearchResultModel documents

### [Materials Electrodes](#/Materials%20Electrodes)

Route providing computed electrode data for a material following the methodology discussed by [Shen _et al._](https://doi.org/10.1038/s41524-020-00422-3)
 Note that electrode data has not been calculated for all materials in the Materials Project database. See the `InsertionElectrodeDoc` and `ConversionElectrodeDoc` schema for a full list of fields returned by this route.

GET

[/materials/insertion\_electrodes/](#/Materials%20Electrodes/search_materials_insertion_electrodes__get)

Get InsertionElectrodeDoc documents

GET

[/materials/conversion\_electrodes/](#/Materials%20Electrodes/search_materials_conversion_electrodes__get)

Get ConversionElectrodeDoc documents

### [Materials Oxidation States](#/Materials%20Oxidation%20States)

Route providing computed oxidation state data for a material following the methodology employed by the [BVAnalyzer](https://pymatgen.org/pymatgen.analysis.bond_valence.html)
 in Pymatgen. Note that oxidation state data has not been calculated for all materials in the Materials Project database. See the `OxidationStateDoc` schema for a full list of fields returned by this route.

GET

[/materials/oxidation\_states/](#/Materials%20Oxidation%20States/search_materials_oxidation_states__get)

Get OxidationStateDoc documents

### [Materials Provenance](#/Materials%20Provenance)

Route providing provenance data for a material such as whether it is theoretical, its associated ICSD entries, and relevant references in literature. Note that provenance data may not be available for all materials in the Materials Project database. See the `ProvenanceDoc` schema for a full list of fields returned by this route.

GET

[/materials/provenance/](#/Materials%20Provenance/search_materials_provenance__get)

Get ProvenanceDoc documents

### [Materials Charge Density](#/Materials%20Charge%20Density)

Route providing computed charge density data for a material following the methodology discussed by [Shen _et al._](https://arxiv.org/abs/2107.03540)
. Please email [heavy.api.use@materialsproject.org](mailto:heavy.api.use@materialsproject.org)
 if you would like to retrieve a large amount of this data. Note that charge densities may not be calculated for all materials in the Materials Project database. See the `ChgcarDataDoc` schema for a full list of fields returned by this route.

GET

[/materials/charge\_density/](#/Materials%20Charge%20Density/search_materials_charge_density__get)

Get ChgcarDataDoc documents

GET

[/materials/charge\_density/{fs\_id}/](#/Materials%20Charge%20Density/get_by_key_materials_charge_density__fs_id___get)

Get a ChgcarDataDoc document by by fs\_id

### [Materials Alloys](#/Materials%20Alloys)

Route for retrevial of information about which hypothetical alloy(s) a given material might belong to, following the methodolgy discussed by [Woods-Robinson, Horton and Persson](https://arxiv.org/pdf/2206.10715)
.

GET

[/materials/alloys/](#/Materials%20Alloys/search_materials_alloys__get)

Get AlloyPairDoc documents

### [Materials Bonds](#/Materials%20Bonds)

Route for retrevial of bonding information for a given material.

GET

[/materials/bonds/](#/Materials%20Bonds/search_materials_bonds__get)

Get BondingDoc documents

### [MPComplete](#/MPComplete)

Route for submitting structures to the Materials Project. If calculations are run with the submitted structure, the submitter will be credited with the submitted public name and email.

GET

[/mpcomplete/](#/MPComplete/search_mpcomplete__get)

Get MPCompleteDoc data

POST

[/mpcomplete/](#/MPComplete/post_data_mpcomplete__post)

Post MPCompleteDoc data

GET

[/mpcomplete/{submission\_id}/](#/MPComplete/get_by_key_mpcomplete__submission_id___get)

Get By Key

### [DOIs](#/DOIs)

Route providing DOI and bibtex reference information for a material. Note that this data may not be available for all materials in the Materials Project database. See the `DOIDoc` schema for a full list of fields returned by this route.

GET

[/doi/](#/DOIs/search_doi__get)

Get DOIDoc documents

### [Molecules Tasks](#/Molecules%20Tasks)

Route for basic task information for DFT calculations in the Materials Project molecules database. Multiple calculations can ultimately be associated with a unique molecule, and are the source of its reported properties. The unique identifier for a calculation is its `task_id`. See the `TaskDocument` schema for a full list of fields returned by this route.

### [Associated Molecules](#/Associated%20Molecules)

Route for 'associated' molecule data. Construction of the Materials Project molecules database occurs in two stages. In the first stage, calculations using the exact same formula, charge, spin multiplicity, and molecular geometry (defined by bond lengths, angles, etc.) are associated. In the second stage, multiple 'associated molecules' with the same basic properties (formula, charge, spin) and connectivity (based on molecular graph isomorphism) are collected, forming the 'core' molecules collection. This route provides access to data for individual 'associated molecules'. The 'Core Molecules' route (/molecules/molecules/) contains data for core molecules. See the `MoleculeDoc` schema for a full list of fields returned by this route.

### [Core Molecules](#/Core%20Molecules)

Route for 'core' molecule data. Construction of the Materials Project molecules database occurs in two stages. In the first stage, calculations using the exact same formula, charge, spin multiplicity, and molecular geometry (defined by bond lengths, angles, etc.) are associated. In the second stage, multiple 'associated molecules' with the same basic properties (formula, charge, spin) and connectivity (based on molecular graph isomorphism) are collected, forming the 'core' molecules collection. This route provides access to data for individual 'associated molecules'. The 'Associated Molecules' route (/molecules/assoc/) contains data for 'associated' molecules. See the `MoleculeDoc` schema for a full list of fields returned by this route.

### [Molecules Partial Charges](#/Molecules%20Partial%20Charges)

Route for molecular partial charge data. See the `PartialChargesDoc` schema for a full list of fields returned by this route.

### [Molecules Partial Spins](#/Molecules%20Partial%20Spins)

Route for molecular partial spin data. See the `PartialSpinsDoc` schema for a full list of fields returned by this route.

### [Molecules Bonds](#/Molecules%20Bonds)

Route for molecular bonding data. See the `MoleculeBondingDoc` schema for a full list of fields returned by this route.

### [Molecules Metal Binding](#/Molecules%20Metal%20Binding)

Route for data regarding metal binding to molecules. See the `MetalBindingDoc` schema for a full list of fields returned by this route.

### [Molecules Orbitals](#/Molecules%20Orbitals)

Route for molecular orbital information obtained via Natural Bonding Orbital analysis. See the `OrbitalDoc` schema for a full list of fields returned by this route.

### [Molecules Redox](#/Molecules%20Redox)

Route for molecular redox information (e.g. ionization energy, reduction free energy, redox potentials). See the `RedoxDoc` schema for a full list of fields returned by this route.

### [Molecules Thermo](#/Molecules%20Thermo)

Route for molecular thermochemistry information. See the `MoleculeThermoDoc` schema for a full list of fields returned by this route.

### [Molecules Vibrations](#/Molecules%20Vibrations)

Route for molecular normal mode and IR spectroscopy data. See the `VibrationDoc` schema for a full list of fields returned by this route.

### [Molecules Summary](#/Molecules%20Summary)

Route for a summary of all data calculated on 'core' molecules in the Materials Project molecules database. See the `MoleculeSummaryDoc` schema for a full list of fields returned by this route.

GET

[/molecules/summary/](#/Molecules%20Summary/search_molecules_summary__get)

Get MoleculeSummaryDoc documents

### [JCESR Electrolyte Genome](#/JCESR%20Electrolyte%20Genome)

Route providing computed data for a legacy molecule such as charge, electron affinity, and ionization energy. The unique identifier for a molecule is its `task_id` (e.g. `mol-45807`). See the `MoleculesDoc` schema for a full list of fields returned by this route.

GET

[/molecules/jcesr/](#/JCESR%20Electrolyte%20Genome/search_molecules_jcesr__get)

Get MoleculesDoc documents

### [Materials Absorption](#/Materials%20Absorption)

GET

[/materials/absorption/](#/Materials%20Absorption/search_materials_absorption__get)

Get AbsorptionDoc documents

### [Materials Chemical Environment](#/Materials%20Chemical%20Environment)

GET

[/materials/chemenv/](#/Materials%20Chemical%20Environment/search_materials_chemenv__get)

Get ChemEnvDoc documents

### [Defect Tasks](#/Defect%20Tasks)

GET

[/defects/tasks/](#/Defect%20Tasks/search_defects_tasks__get)

Get DefectTaskDoc documents

#### Schemas

AbsorptionDoc

Expand all**object**

AlloyPairDoc

Expand all**object**

AnalysisDoc

Expand all**object**

Author

Expand all**object**

BSObjectDoc

Expand all**object**

BSPathType

Expand all**string**

BandStructureSummaryData

Expand all**object**

BandstructureData

Expand all**object**

BatteryType

Expand all**string**

BlessedCalcs

Expand all**object**

BondingComposite

Expand all**object**

BondingDoc

Expand all**object**

BulkModulus

Expand all**object**

Calculation

Expand all**object**

CalculationInput

Expand all**object**

CalculationOutput

Expand all**object**

ChemEnvDoc

Expand all**object**

ChgcarDataDoc

Expand all**object**

ComplianceTensorDoc

Expand all**object**

Component

Expand all**object**

CondensedStructureData

Expand all**object**

Conditions

Expand all**object**

ConversionElectrodeDoc

Expand all**object**

ConversionVoltagePairDoc

Expand all**object**

CrystalSystem

Expand all**string**

CustodianDoc

Expand all**object**

DOIDoc

Expand all**object**

DOSObjectDoc

Expand all**object**

DOSProjectionType

Expand all**string**

Database

Expand all**string**

DecompositionProduct

Expand all**object**

DefectTaskDoc

Expand all**object**

DeprecationDoc

Expand all**object**

DeprecationMessage

Expand all**string**

DielectricDoc

Expand all**object**

DosData

Expand all**object**

DosSummaryData

Expand all**object**

EOSDoc

Expand all**object**

Edge

Expand all**string**

ElasticTensorDoc

Expand all**object**

ElasticityDoc

Expand all**object**

ElectronPhononDisplacedStructures

Expand all**object**

ElectronicStep

Expand all**object**

ElectronicStructureDoc

Expand all**object**

Element

Expand all**string**

EmmetMeta

Expand all**object**

EntriesCompositionSummary

Expand all**object**

EntryDoc

Expand all**object**

Error

Expand all**object**

ExtractedMaterial

Expand all**object**

FermiDoc

Expand all**object**

FindStructure

Expand all**object**

FittingData

Expand all**object**

FormulaAutocomplete

Expand all**object**

FormulaPart

Expand all**object**

FrequencyDependentDielectric

Expand all**object**

GBSearchData

Expand all**object**

GBTypeEnum

Expand all**string**

GrainBoundaryDoc

Expand all**object**

HTTPValidationError

Expand all**object**

History

Expand all**object**

InputDoc

Expand all**object**

InsertionElectrodeDoc

Expand all**object**

InsertionVoltagePairDoc

Expand all**object**

IonicStep

Expand all**object**

LevelOfTheory

Expand all**string**

MPCompleteDataStatus

Expand all**string**

MPCompleteDoc

Expand all**object**

MagnetismDoc

Expand all**object**

MaterialsDoc

Expand all**object**

Meta

Expand all**object**

MetalBindingComposite

Expand all**object**

MetalBindingData

Expand all**object**

MineralData

Expand all**object**

MoleculeSummaryDoc

Expand all**object**

MoleculesDoc

Expand all**object**

MultipolesComposite

Expand all**object**

Operation

Expand all**object**

OperationTypeEnum

Expand all**string**

OrbitalComposite

Expand all**object**

OrbitalType

Expand all**integer**

Ordering

Expand all**string**

OrigInputs

Expand all**object**

OutputDoc

Expand all**object**

OxidationStateDoc

Expand all**object**

PartialChargesComposite

Expand all**object**

PartialSpinsComposite

Expand all**object**

PhaseDiagramDoc

Expand all**object**

PhononBSDOSDoc

Expand all**object**

PiezoelectricDoc

Expand all**object**

PointGroupData

Expand all**object**

Potcar

Expand all**object**

PotcarSpec

Expand all**object**

PropertyOrigin

Expand all**object**

ProvenanceDoc

Expand all**object**

ReactionFormula

Expand all**object**

RedoxComposite

Expand all**object**

Response\[AbsorptionDoc\]

Expand all**object**

Response\[AlloyPairDoc\]

Expand all**object**

Response\[BSObjectDoc\]

Expand all**object**

Response\[BondingDoc\]

Expand all**object**

Response\[ChemEnvDoc\]

Expand all**object**

Response\[ChgcarDataDoc\]

Expand all**object**

Response\[ConversionElectrodeDoc\]

Expand all**object**

Response\[DOIDoc\]

Expand all**object**

Response\[DOSObjectDoc\]

Expand all**object**

Response\[DefectTaskDoc\]

Expand all**object**

Response\[DeprecationDoc\]

Expand all**object**

Response\[DielectricDoc\]

Expand all**object**

Response\[EOSDoc\]

Expand all**object**

Response\[ElasticityDoc\]

Expand all**object**

Response\[ElectronicStructureDoc\]

Expand all**object**

Response\[EntryDoc\]

Expand all**object**

Response\[FermiDoc\]

Expand all**object**

Response\[FindStructure\]

Expand all**object**

Response\[FormulaAutocomplete\]

Expand all**object**

Response\[GrainBoundaryDoc\]

Expand all**object**

Response\[InsertionElectrodeDoc\]

Expand all**object**

Response\[MPCompleteDoc\]

Expand all**object**

Response\[MagnetismDoc\]

Expand all**object**

Response\[MaterialsDoc\]

Expand all**object**

Response\[MoleculeSummaryDoc\]

Expand all**object**

Response\[MoleculesDoc\]

Expand all**object**

Response\[OxidationStateDoc\]

Expand all**object**

Response\[PhaseDiagramDoc\]

Expand all**object**

Response\[PhononBSDOSDoc\]

Expand all**object**

Response\[PiezoelectricDoc\]

Expand all**object**

Response\[ProvenanceDoc\]

Expand all**object**

Response\[RobocrystallogapherDoc\]

Expand all**object**

Response\[SimilarityDoc\]

Expand all**object**

Response\[SubstratesDoc\]

Expand all**object**

Response\[SummaryDoc\]

Expand all**object**

Response\[SummaryStats\]

Expand all**object**

Response\[SurfacePropDoc\]

Expand all**object**

Response\[SynthesisSearchResultModel\]

Expand all**object**

Response\[TaskDoc\]

Expand all**object**

Response\[ThermoDoc\]

Expand all**object**

Response\[TrajectoryDoc\]

Expand all**object**

Response\[XASDoc\]

Expand all**object**

RobocrystallogapherDoc

Expand all**object**

RunStatistics

Expand all**object**

RunType

Expand all**string**

ShearModulus

Expand all**object**

SimilarityDoc

Expand all**object**

SimilarityEntry

Expand all**object**

SoundVelocity

Expand all**object**

Spin

Expand all**integer**

Status

Expand all**string**

SubstratesDoc

Expand all**object**

SummaryDoc

Expand all**object**

SummaryStats

Expand all**object**

SurfaceEntry

Expand all**object**

SurfacePropDoc

Expand all**object**

SymmetryData

Expand all**object**

SynthesisSearchResultModel

Expand all**object**

SynthesisTypeEnum

Expand all**string**

TaskDoc

Expand all**object**

TaskState

Expand all**string**

ThermalConductivity

Expand all**object**

ThermoComposite

Expand all**object**

ThermoDoc

Expand all**object**

ThermoType

Expand all**string**

TrajectoryDoc

Expand all**object**

Type

Expand all**string**

ValidationError

Expand all**object**

Value

Expand all**object**

Values

Expand all**object**

VaspObject

Expand all**string**

VibrationComposite

Expand all**object**

XASDoc

Expand all**object**

XASSearchData

Expand all**object**

CalcType

Expand all**string**

TaskType

Expand all**string**

CalcType

Expand all**string**

TaskType

Expand all**string**