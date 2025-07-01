# https://docs.materialsproject.org

This is public documentation for the [Materials Project](https://materialsproject.org) (MP). The Materials Project is a decade-long effort from the Department of Energy to pre-compute properties of "materials" and make this data publicly available, with the intent of accelerating the process of materials discovery. In this context, a material can mean either an inorganic crystal (like silicon), or a molecule (like ethylene). Possible applications are vast, but might include better batteries, solar energy, water splitting, optoelectronics, catalysts and more (see [here](https://materialsproject.org/about/publications) for a list of publications).

## [Direct link to heading](\#table-of-contents)    Table of Contents

### [Direct link to heading](\#methodology)    Methodology

This section contains information on how we generate and validate our computed data sets.

[Methodology](/methodology/materials-methodology)

### [Direct link to heading](\#apps)    Apps

This section talks about how we present information on the Materials Project website as "apps", and what these apps contain.

### [Direct link to heading](\#getting-involved)    Getting Involved

The Materials Project is a public, collaborative project, offered free of charge. It only succeeds thanks to the efforts of everyone who participates! Visit this section to learn how to get involved.

### [Direct link to heading](\#errors)    Errors

If you notice an error or omission, please let us know at [our user forum located at matsci.org/materials-project](https://matsci.org/materials-project).

The Materials Project documentation is a living document and always a work in progress.

[NextFrequently Asked Questions (FAQ)](/frequently-asked-questions)

Last updated 1 month ago

---

# https://docs.materialsproject.org/frequently-asked-questions

This page contains answers to common questions about the Materials Project.

See also our "Glossary of Terms" page which defines common terms in use by Materials Project.

[Glossary of Terms](/frequently-asked-questions/glossary-of-terms)

## [Direct link to heading](\#how-do-i-sign-in-to-the-materials-project)    How do I sign in to the Materials Project?

You can login to the Materials Project either using an existing social identity provider (currently GitHub, Google, Facebook, Microsoft or Amazon) or via an email link.

Be aware, your Materials Project account is linked to both your email address and the method that you log in. If you log in via a different method, this will be registered as a new account.

Here are some issues people have encountered when trying to sign in the [Materials Project](https://materialsproject.org/) website, and their solutions:

- **I want to log in with my social identity provider (GitHub/Google/Facebook/Microsoft/Amazon), but I can’t.**



Ensure that your password for your provider is correct (go to their site and log in there), ensure that you have a full name set on that account, and ensure that you allow Materials Project to see your basic profile info (name and email address).



You also may be behind a firewall that doesn’t allow GitHub/Google/Facebook/Microsoft/Amazon. In that case, use our email based option instead.

- **I appear to sign in OK (the popup goes away), but then I remain on the sign-in screen.**



It may take a few seconds, depending on your connection, to actually get logged in. This is because we have an external identity provider verify your email address so that we don’t have to store any passwords on our servers.



You may also have an older browser that won’t work well with our website at all. The latest version of Mozilla Firefox, Google Chrome, or Microsoft Edge will work well. Older versions of Internet Explorer will not work.

- **I tried using the email option several times but haven’t received a login** link.



We currently don’t do any validation of your email addresses, so if it “looks right”, i.e. you mistype [_myname@gmail.com_](mailto:myname@gmail.com) as _myname@_ _**gmali**_ _.com_, we will still try to send to the wrong address. Also check your "Spam" or "Junk" folder in case the login email has been flagged.



There is a known issue with Tencent @qq.com addresses, where Tencent throttles delivery and you might not get an email within a reasonable amount of time. Please consider using an alternative to your @qq.com address for login.


## [Direct link to heading](\#how-do-i-cite-materials-project)    How do I cite Materials Project?

Citations are appropriate wherever Materials Project data, methods or output are used. See this page on the Materials Project website for more information:

[![Logo](https://materialsproject.org/assets/favicon.ico?m=1657659155.0)Materials Project - How to CiteMaterials Project](https://materialsproject.org/about/cite)

How to Cite Materials Project

There is a canonical Materials Project citation, and additional citations for specific properties or tools. See also the [Database Versions](/changes/database-versions) page for information on how to cite a specific database version.

## [Direct link to heading](\#where-do-the-material-properties-shown-on-materials-project-come-from)    Where do the material properties shown on Materials Project come from?

The Materials Project core data is all calculated in-house by the Materials Project team using a variety of simulation methods. To understand the quality of these predictions, it is crucial to read the peer-reviewed publications from the Materials Project where each property is benchmarked as much as possible against known experimental values: this will give an estimate of typical error and, importantly, any systematic error that may be present.

## [Direct link to heading](\#why-are-the-lattice-parameters-different-to-what-i-expect)    Why are the lattice parameters different to what I expect?

The same crystal structure can have multiple, equivalent sets of lattice parameters depending on what crystallographic "setting" is used.

Typically, there are two sets of lattice parameters reported. Lattice parameters can be defined for the **primitive** **cell**, which is a definition of the crystal with the fewest number of atoms and therefore convenient for simulations and other uses, and the **conventional cell**, which is typically easier to visualize and more like you will see in textbooks.

If the lattice parameters are very different to what you expect, check the setting first!

Some systematic errors are also present. These will typically be an over-estimation of 1–3% for most crystals. Layered crystals will also typically have significant error in the interlayer distances since van der Waals interactions are not well-described by the simulation methods (PBE) used by Materials Project. These systematic errors will be improved as Materials Project switches to user newer simulation methods (r2SCAN). See [Calculation Details](/methodology/materials-methodology/calculation-details) for more information.

## [Direct link to heading](\#why-is-the-band-gap-different-to-what-i-expect)    Why is the band gap different to what I expect?

Electronic band gaps are difficult to calculate reliably from first principles, especially using methods that scale well to hundreds of thousands of materials. The method used by the Materials Project (PBE) _systematically underestimates_ band gaps.

While it would be possible to provide higher quality calculations for a select number of materials, with more accurate band gaps, it is noted that for materials discovery purposes it is useful to have a dataset that has the same systematic error. See Electronic Structure for more information.

## [Direct link to heading](\#why-has-a-value-changed-on-materials-project)    Why has a value changed on Materials Project?

The Materials Project presents the data it generates in two ways:

1. As individual calculations. These are always the same, and as far as possible Materials Project tries to ensure all historical calculations remain available. Typically, only advanced users will access information about individual calculations.

2. As aggregated information. This is information generated from a combination of individual calculations. This information is what is presented on the public "material details" pages, and is what most users will access. As new, improved calculations are performed, this aggregated information can change.


The Materials Project periodically updates this aggregated information in the form of new database releases. See [Database Versions](/changes/database-versions) for information on the latest database releases.

If performing scientific research with Materials Project data, make sure to cite the database version from which the data was retrieved. See [How to Cite](https://materialsproject.org/about/cite) for more information.

## [Direct link to heading](\#what-is-a-task_id-and-what-is-a-material_id-and-how-do-they-differ)    What is a "task\_id" and what is a "material\_id" and how do they differ?

Every database needs a unique key which can be used to distinguish one entry from another. In the Materials Project, each unique material is given a `material_id` (also referred to in various places as mp-id, mpid, MPID). This allows a specific polymorph of a given material to be referenced. For example, wurtzite GaN is assigned the `material_id` of [`mp-804`](https://materialsproject.org/materials/mp-804), while zinc blende GaN is assigned a `material_id` of [`mp-830`](https://materialsproject.org/materials/mp-830).

### [Direct link to heading](\#how-does-a-material_id-get-assigned)    How does a "material\_id" get assigned?

The Materials Project is a computational resource. All of the information on a given _material details page_ is actually a combination of data generated from many individual calculations or "tasks". It is also important that these tasks also have unique identifiers.

When a task is added to the Materials Project database, it will get an identifier assigned with the format `mp-[0-9]` ("mp-" with numbers after it). These identifiers are assigned sequentially, so smaller numbers usually refer to older calculations. An identifier referring to an _individual calculation task_ are known as a `task_id`.

When the Materials Project database is built, a unique material will then have a collection of multiple different task\_ids associated with it. The numerically smallest `task_id` will then become the `material_id`. This ensures that, as new, additional calculations are associated with the same material, its `material_id` should not change.

### [Direct link to heading](\#in-the-past-i-have-seen-material_ids-that-start-with-mvc-what-are-these)    In the past, I have seen material\_ids that start with "mvc", what are these?

Some calculation tasks were associated with a search for multivalent cathode materials. These tasks were given the prefix `mvc-` instead of `mp-` and thus some materials also had the prefix `mvc-`. However, this caused confusion and this approach has been retired. Tasks with the prefix `mvc-` still exist since the `task_id` cannot change, but a `material_id` will now always start with an `mp-` prefix by convention provided that at least one task associated with that material has the `mp-` prefix.

### [Direct link to heading](\#do-material_ids-ever-change-do-task_ids-ever-change)    Do material\_ids ever change? Do task\_ids ever change?

A `task_id` will never change. It will always refer to the same, individual calculation task.

A `material_id` might change in rare instances, such as the removal of the `mvc-` prefix, although this is avoided wherever possible.

If a `material_id` _does_ change, we ensure a redirect on the website is always in place, and the new `material_id` can also be found programmatically with the API using the `get_material_id_from_task_id()` function. This way, any publications or research that reference an older `material_id` are still valid, and the relevant data can still be retrieved.

## [Direct link to heading](\#what-does-______-mean)    What does \_\_\_\_\_\_ mean?

Consult our glossary here:

[Glossary of Terms](/frequently-asked-questions/glossary-of-terms)

If a term is used in Materials Project but is not listed, [let us know](/community/getting-help) and we will add it.

[PreviousIntroduction](/) [NextGlossary of Terms](/frequently-asked-questions/glossary-of-terms)

Last updated 11 months ago

---

# https://docs.materialsproject.org/frequently-asked-questions/glossary-of-terms

**Builder.** A builder is a little script written in the Python programming language that helps create new database collection(s) from input database collection(s). It's typically used to allow common analysis tasks to be repeated automatically, for example the calculation of "energies above hull" when new calculations are added to the database. Builders are an essential step in the Materials Project database release process and are formalized with the `emmet` code.

**Chemical system.** On Materials Project, a chemical system is a set of materials whose members all contain the same elements. It is usually noted with as dash-delimited list of elements. For example, the "Ga-In-N" chemical system would contain all materials containing Ga, In or N or combinations of these elements (Ga, In, N2, GaN, InGaN, etc.).

**Correction scheme.** The Materials Project performs calculations using a simulation technique with known systematic errors. A correction scheme is employed to adjust energies based on the elements present in a material to address these systematic errors. Only elements for which sufficient experimental data is available can be corrected.

**Energy above hull.** A measure of a material's thermodynamic stability. This value refers to a mathematical construction that can be calculated from a set of formation energies and compositions known as a convex hull, and often referred to here as a "phase diagram." However, unlike most phase diagrams, convex hulls are usually given without a temperature axis since the simulation technique used (DFT) gives predictions at zero temperature. A material which lies "on the convex hull" is predicted to be thermodynamically stable, while off the hull is predicted to be metastable or unstable. Values above 200 meV/atom are considered very large and suggest an unstable material that might not be synthesizable, however this ceiling differs significantly by chemistry. Energies above hull are given as a guide and subject to both limits of calculation precision (several meV) and also of calculation accuracy due to limitations of the simulation technique used, where errors can be significant in certain chemistries.

**Mixing scheme.** The Materials Project uses two slightly different simulation techniques depending on the elements present in a material. These are GGA (Generalized Gradient Approximation) and GGA+U, where the +U (Hubbard correction) is a correction applied to address systematic deficiencies in GGA when simulating elements with highly localized electrons such as d-orbitals or f-orbitals. Energies from these respective techniques are not directly comparable with each other, so a mixing scheme is employed such that elements can be compared. Details of the mixing scheme can be found in [this paper](https://doi.org/10.1103/PhysRevB.84.045115).

[PreviousFrequently Asked Questions (FAQ)](/frequently-asked-questions) [NextChanges and Updates](/changes)

Last updated 11 months ago

---

# https://docs.materialsproject.org/changes

The Materials Project is an active, academic research project. Changes are common as new research methods become available, and the quality and kind of data we present changes, and also as a result of organizational needs. This page summarizes major changes in different aspects of the Materials Project.

## [Direct link to heading](\#upcoming-changes)    Upcoming Changes

This documentation will continue to be improved. New documentation is currently being written for each of the Materials Project "apps". Some pages may be blank until this is completed.

## [Direct link to heading](\#previous-changes)    Previous Changes

### [Direct link to heading](\#database)    Database

The Materials Project database is constantly evolving as new and better calculations become available, both as a result of new features and better methods, and also as errors or problems are identified and fixed.

See the following documentation page for a list of changes to the Materials Project database:

[Database Versions](/changes/database-versions)

### [Direct link to heading](\#api)    API

The Materials Project API has recently undergone a significant modernization effort. The new Materials Project website is exclusively powered by this API.

See the following documentation page for more information:

[Differences between new and legacy API](/downloading-data/differences-between-new-and-legacy-api)

### [Direct link to heading](\#website)    Website

The Materials Project has recently undergone a major change in its website architecture. More information on this can be seen in the [release announcement](https://medium.com/materials-project/announcing-a-new-materials-project-2628ded751c).

It is recommended that the URL [https://materialsproject.org](https://materialsproject.org) is used as the primary location of the Materials Project website, however a specific website version can be visited via the following links:

- [https://next-gen.materialsproject.org](https://next-gen.materialsproject.org) will always take visitors to the latest Materials Project website with the newest database version available.

- [https://legacy.materialsproject.org](https://legacy.materialsproject.org) will take visitors to a frozen snapshot of the older Materials Project website. This is powered by an older version of the database with known issues. The legacy website is being left online for some time as we fully transition to the next-gen website, and to allow users time to make any adjustments as necessary for features that may only be available on the legacy website, however the legacy website will be taken offline in due course.


See the website changelog for a detailed list of recent changes:

[Website Changelog](/changes/website-changelog)

### [Direct link to heading](\#documentation)    Documentation

The Materials Project documentation has gone through several iterations, powered previously by MediaWiki and MkDocs software. The current version is powered by GitBook. This switch was made to allow more easy and rapid changes to the documentation, in the hopes of ensuring documentation is maintained at a consistent, high quality.

The current documentation is also available via GitHub at [https://github.com/materialsproject/public-docs](https://github.com/materialsproject/public-docs). Edits and improvements from external users are very welcome, please submit a "pull request" with any suggest change or use the "Edit in GitHub" button on the relevant page.

The previous MkDocs documentation is [still available](https://github.com/materialsproject/docs) for the historical record, and the older MediaWiki documentation are currently offline but available on request. However, the current version of the documentation should contain all necessary information including historical information. An effort has been made to ensure URLs remain the same during the transition from the previous MkDocs-powered documentation to the new GitBook-powered documentation.

[PreviousGlossary of Terms](/frequently-asked-questions/glossary-of-terms) [NextDatabase Versions](/changes/database-versions)

Last updated 2 years ago

---

# https://docs.materialsproject.org/changes/database-versions

This page contains a summary of major changes for each version of the Materials Project database.

We are aware of a community need for more detailed change logs, and hope to improve our reporting for future database versions.

Database versions are labelled via the date they become generally available to the public. For advanced users, the public database label is mapped to an internal pull request number, and this is accessible in the API via the `builder_meta` key.

You can verify the current database version powering the website on the footer of every page. If you are using the API, there is a `get_database_version` method available.

## [Direct link to heading](\#v2023.11.1)    v2023.11.1

- Improved and expanded set of elasticity data. Note that there are schema changes with how it is accessed in `SummaryDoc` and `ElasticityDoc`.

- Conversion electrode data added alongside existing insertion electrode data.

- ~10k new materials added, with ~5k deprecated. This includes a temporary deprecation of all compounds containing `Yb` while they are being re-run. This is in response to pseudopotential issues identified which were providing incorrect energies.


## [Direct link to heading](\#v2022.10.28)    v2022.10.28

This database build incorporates Materials Project’s (R2)SCAN calculations as pre-release data. The default fields returned by the website and API will remain unchanged from the previous release at the GGA(+U) level of theory, but the (R2)SCAN data is now available for advanced users. Either see the “Pre-release Data” section of a relevant material details page, generate an R2(SCAN) phase diagram with the Phase Diagram app, or access the data via the thermo API endpoint. This database release also incorporates several new perovskite materials from a collaboration with Zachary Bare, University of Colorado.

## [Direct link to heading](\#v2021.11.10)    v2021.11.10

This will be the first release with our new website and API. It does not contain any new data but is built using our new database building methods and is largely consistent with the previous database release. Some changes exist to the previous release due to improvements to detection of multi-anion systems leading to changes in the applied formation energy corrections.

Be aware, database version v2021.11.10 onwards is _only_ available on the new Materials Project website and API. The [legacy website](https://legacy.materialsproject.org) and [legacy API](/downloading-data/differences-between-new-and-legacy-api) are frozen to the v2021.05.13 database release.

## [Direct link to heading](\#v2021.05.13)    v2021.05.13

This release updates the energy correction scheme we use to generate phase diagrams and compute formation energies. **As with any new database release, formation energies for many compounds have changed**; however in this case the change is due only to our new energy correction scheme and not to any new data. We are proud to report that the new correction scheme has reduced the overall error in formation energy in our database by 7% compared to experiment.

You can see details of each correction that has been applied by inspecting the `energy_adjustments` attribute of a `ComputedEntry` retrieved via the API. In addition, the new correction scheme is available for manual use via the `MaterialsProject2020Compatibility` class in pymatgen.

We realize that this change may be disruptive to ongoing work, and want to assure you that the historical corrections are still available in pymatgen if needed. They may be recovered by manually reprocessing `ComputedEntry` using the legacy `MaterialsProjectCompatibility` class. An example notebook demonstrating how to do this available [on matgenb 25](https://github.com/materialsvirtuallab/matgenb/blob/3dc1e275f0a9ceadc032d83d71601676530d736e/notebooks/2021-5-12-Explanation%20of%20Corrections.ipynb).

Below we summarize the most significant changes associated with the new `MaterialsProject2020Compatibility` correction scheme. For complete details and documentation, please refer to [this manuscript 32](https://chemrxiv.org/articles/preprint/A_Framework_for_Quantifying_Uncertainty_in_DFT_Energy_Corrections/14593476).

**1\. Refitted corrections for legacy species**
Corrections applied to oxygen compounds, diatomic gases, and transition metal oxides and fluorides have been refit using more up to date DFT calculations and a larger compilation of computed and experimental formation enthalpy data.

**2\. Corrections for additional species**
We have added corrections for Br, I, Se, Si, Sb, and Te, which did not previously have energy corrections. As a result, formation energies for materials containing these species will generally be lower than they were previously.

**3\. Diatomic gas corrections moved to compounds**
Previously, corrections for H, F, Cl, and N were applied to the elements. One consequence of this was that polymorphs of H2, N2, Cl2 and F2 were _always_ assigned a zero energy above hull, even if some polymorphs were higher in energy. This made interpretation of these values confusing. With this release, energy corrections are applied to the material (e.g., LiH) and not the element. This also means that unstable polymorphs of diatomic gases will now have non-zero `e_above_hull`

**4\. Oxidation state based corrections**
Our build process now estimates the likely oxidation states of each species in a material, and uses this information to intelligently apply corrections to anionic species only when their estimated oxidation state is negative. For example, in the compound `MoCl3O`, estimated oxidation states for both Cl and O are negative, so both anions receive corrections.

Our algorithms are not always successful in predicting the oxidation state. When this occurs, we apply anion corrections to only the most electronegative element in the material. As a result, some ternary or higher compounds in the database may be destabilized in this release because their oxidation states could not be determined. This is the case for MoCl5O (mp-1196724) for example, which does not receive a Cl correction because O is more electronegative.

If this affects your work, you can manually assign oxidation states by populating the `oxidation_states` key of the `.data` attribute of any `ComputedEntry` and then reprocessing the data using `MaterialsProject2020Compatibility`.

**5\. Uncertainty Quantification**
We now compute the estimated uncertainty associated with the energy corrections on a material. Uncertainties reflect the measured uncertainty in the underlying experimental data that we use to determine the corrections, as well as uncertainty associated with the fitting procedure itself. This information enables new methods of assessing phase stability, as described in [this manuscript 32](https://chemrxiv.org/articles/preprint/A_Framework_for_Quantifying_Uncertainty_in_DFT_Energy_Corrections/14593476)

**A Note for API and MPRester Users**

For API users, if you are retrieving formation energies directly via the API, you will get the correct, latest formation energies from the current database release. However, if you are using `get_entries` or `get_pourbaix_entries` which apply the correction scheme on-the-fly, make sure to update to the latest version of pymatgen (v2022.0.8 or later) to get the correct values. If you are using pymatgen v2021 or earlier, this will use the old correction scheme by default when using `get_entries` and `get_pourbaix_entries`.

## [Direct link to heading](\#v2021.03.22)    v2021.03.22

This release updates some older materials with new calculations, and adjusts our rules for deprecating older calculations. It does not contain any new materials. Thanks to the new calculations many materials that were previously deprecated are now accessible again. This release is in preparation for a switch to our new compatibility scheme which will improve our predictions of formation energy.

## [Direct link to heading](\#v2021.02.08)    v2021.02.08

We had a small new database release today, this introduces new higher-quality calculations for around 30,000 materials. It also deprecates 78 materials since we currently do not have calculations for these materials that match our current quality standards; we hope to restore these 78 materials in a subsequent release. For an exact list, please see the attached file.

As a reminder, all historical calculation tasks remain available via our API and the task detail pages, and information on deprecated materials also remain available via the API. More information on our deprecation policy is in our documentation. We continue our work on better ways communicate database diffs and to more easily provide access to historical information, so stay tuned for future announcements here.

[db\_v2020\_09\_08\_to\_v2021\_02\_08\_diff.yaml](https://matsci.org/uploads/short-url/CITebuwCI46kWUbfjIWfMFf6ye.yaml) (376.9 KB)

## [Direct link to heading](\#v2020.09.08)    v2020.09.08

This releases addresses issues noticed in the previous release with formation energies and updates the energies of approximately 6k materials where this error was greatest. We are planning a further supplemental release.

We’re also looking at ways to put in place a process to be more transparent with database changes and updates to share more specifically what has changed, as well as providing means to access historical versions of the database, since we know this is a common requirement.

Note that, wherever possible, we continue to keep individual historical calculation data available via its _task\_id_ even in cases where the aggregated information (such as that presented on the materials detail page) might change.

**V2020.08.20 Released**

In this release we have added thousands of new band structure and density of states calculations, improving our overall material coverage and data quality. Additionally, we have overhauled the plotting for these quantities on the material details page. This is a first step in improving the electronic structure data within the Materials Project as part of our [new tool set 46](https://www.nature.com/articles/s41524-020-00383-7) for band structure calculations.

We are also working through an on-going issue affecting the energies of a small number of materials. In the previous release, we added a large batch of higher-quality calculations for our energetics as well as fixing numerous bugs. However, we discovered an error in our calculation parameters leading to larger energies than expected for a minority of materials and issues such as those discussed [here 27](https://matsci.org/t/inconsistent-energy-between-mp-756366-and-mp-763752-li/4648). We are currently re-running these calculations and will be fixing this data in a supplemental update in the next few weeks. We advise anybody performing large screening studies to do so with caution or wait until this supplemental update has been released.

## [Direct link to heading](\#v2020.06)    v2020.06

In this database release, we have added several thousand materials and many magnetic ground states, improved the quality of our energetics, and fixed many bugs. This database release is part of on-going efforts in 2020 to improve database reliability and quality, following the introduction of our deprecation process last year. There are still known issues with this release which we are working to address, please let us know if you encounter any in our forum.

## [Direct link to heading](\#v2019.12.05)    **v2019.12.05**

The issue mentioned in 2019-12-04 has now been addressed, however approximately 7% of materials saw errors in their reported energies above hull of greater than 0.05 eV/atom. Values calculated via _pymatgen_ or via the phase diagram app on the website during this time were correct, while values reported on the materials details page and via the `e_above_hull` API key were incorrect.

We encourage users who accessed convex hulls from the website between the latest database release and 2019-12-05 to re-check any values obtained from the website.

We apologize for the error, and will be incorporating additional checks into our automated testing to prevent similar errors in the future.

## [Direct link to heading](\#v2019.12.04)    **v2019.12.04**

We are aware of an on-going issue with the reported energies above hull on the materials detail pages. We will update this thread when a fix has been fully implemented and with further details.

Until this issue is fully resolved, correct energies above hull can be retrieved using _pymatgen_ as follows:

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
from pymatgen import MPRester
from pymatgen.analysis.phase_diagram import PhaseDiagram

with MPRester(YOUR_API_KEY_HERE) as mpr:
    # replace with your elements and mp-id of interest
    entries = mpr.get_entries_in_chemsys(['Li','Co', 'O'])
    entry_of_interest = mpr.get_entry_by_material_id('mp-19128')

phase_diagram = PhaseDiagram(entries)
e_above_hull = phase_diagram.get_e_above_hull(entry_of_interest)

print("e_above_hull", e_above_hull)
```

If you have not previously used _pymatgen_, it is a Python code and can be installed using `pip install pymatgen` or `conda install --channel conda-forge pymatgen`.

_Note: the above information for v2019-12-04 is now out of date._

## [Direct link to heading](\#v2019.11.21)    v2019.11.21

During deployment of the new v2019.11 database, there was temporary issue with generating interactive phase diagrams leading to incorrect formation enthalpies for a small number of chemical systems. This has now been fixed. Data presented on the materials detail pages was unaffected by this issue.

## [Direct link to heading](\#v2019.11)    v2019.11

- Introduced 3,971 new materials

- Amorphous materials added with `amorphous` tag

- Added `theoretical` which is True when the material matches no known experimental structure from ICSD

- Fixed several inconsistency bugs for `band_gap`, piezo tensors, elastic warnings, and total magnetic moment.


## [Direct link to heading](\#v2019.05)    v2019.05

- Introduced a new `deprecated` field to materials. By default the website and API only search for materials that are not deprecated: {“deprecated”: false}.

- Deprecated 15,000 and added 3,600 new materials. We will be recomputing the deprecated materials to fill these spaces back up. Some of these new relaxations may end up matching current materials, so the total number of materials is not guaranteed to be the same as in V2019.02. This also affects downstream properties. Most notably, ~3k elastic tensors associated with the deprecated materials have been removed from the database and are no longer accessible.

- Fixed an issue with sandboxes not properly building the whole hull. Previously, only the sandboxed chemical systems were being recalculated for energy\_above\_hull searches


## [Direct link to heading](\#v2019.02)    v2019.02

- Added over 47,000 new materials from orderings of disordered ICSD as well as compounds from the Pauling File

- Finalized enforcing symmetry on piezo tensors

- Moved third order elastic data to elasticity\_third\_order so that people are not swamped by the mountain of information associated with it.


## [Direct link to heading](\#v2018.12)    v2018.12

- Adjusted the mp-id naming scheme to fix “mvc” ids taking over old mp-ids.

- Fixed piezoeletric max\_direction to be a miller index rather than a unit vector.


## [Direct link to heading](\#v2018.11)    v2018.11

- Changed the grouping of magnetic materials to aggregate all magnetic orderings of a given material into a single material-id, and report the lowest energy ordering

- Fixed incorrect calculation and display of polycrystalline dielectric constants

- Fixed labeling of all materials as high-pressure. Note we’re parsing ICSD tags for this labeling so while some materials may not conventionally be considered high-pressure, a single matching ICSD entry can tag a material as such. We would love to hear comments on how we could better tag high-pressure materials

- Begun enforcing the symmetry of the structure on piezo tensors. In general, this reduces the expected piezo value.


[PreviousChanges and Updates](/changes) [NextWebsite Changelog](/changes/website-changelog)

Last updated 10 months ago

---

# https://docs.materialsproject.org/changes/website-changelog

#### [Direct link to heading](\#id-2022-12-16-7ca3bcd3)    2022-12-16 (7ca3bcd3)

- Fix issue with API query, see [here](https://matsci.org/t/rest-query-returned-with-error-status-code-500/45793).

- Better integration between `MPRester` and `MPContribs` API python clients.



- **Users of the new API should upgrade to** `mp-api>=0.30.5` **and** `mpcontribs-client>=5.0.4`


#### [Direct link to heading](\#id-2022-12-02-13f229ed)    2022-12-02 (13f229ed)

- Fix an incorrect unit label for elasticity data on the new website. Thank you to Serge Maalouf for reporting.



- Data returned from the API was correct and unaffected by this error.


- Fix for insufficient precision in reporting atomic co-ordinates of some materials. Kindly reported by Branton Campbell for the entry `mp-1106336`.



- Data returned from the API was correct and unaffected by this error.


- An issue with displaying "task detail" pages is resolved.


#### [Direct link to heading](\#id-2022-08-09-f2aa3e0a)    2022-08-09 (f2aa3e0a)

- Resolved a bug with "MOF Explorer" detail pages not loading.

- We are investigating an issue with the "Crystal Toolkit" app.



- This was resolved.


#### [Direct link to heading](\#id-2022-07-28-e7527896)    2022-07-28 (e7527896)

- Added "Alloy Systems" section to the material details pages.



- This is a preview of a new feature and is not yet peer-reviewed.

- More information on the methodology is available [here](https://arxiv.org/abs/2206.10715).

- Examples of this feature might be seen on the materials detail page for [CdTe](https://materialsproject.org/materials/mp-406) or [GaN](https://materialsproject.org/materials/mp-804).


#### [Direct link to heading](\#id-2022-07-12-5d802243)    2022-07-12 (5d802243)

- Fixed an issue with permuted axis labels in the Equations of State plots, kindly reported by [zzyfor2019](https://matsci.org/u/zzyfor2019) on the forum



- The data returned by the API was correct and unaffected by this error


- Fixed an issue with swapped labels in the Battery Explorer, kindly reported by 施荣鑫 via email



- The data returned by the API was correct and unaffected by this error


[PreviousDatabase Versions](/changes/database-versions) [NextDocumentation Credit](/credit)

Last updated 1 year ago

---

# https://docs.materialsproject.org/credit

The Materials Project documentation is a collaborative effort between Materials Project staff, contributors, and researchers including graduate students, postdocs and members of the Materials Project community.

A recent list of contributors can be found here:

[![Logo](https://github.com/fluidicon.png)Contributors to materialsproject/public-docsGitHub](https://github.com/materialsproject/public-docs/graphs/contributors)

Recent list of Materials Project documentation contributors.

See also the "Documentation Authors" sections on individual documentation pages.

[PreviousWebsite Changelog](/changes/website-changelog) [NextGetting Help](/community/getting-help)

Last updated 11 months ago

---

# https://docs.materialsproject.org/community/getting-help

## [Direct link to heading](\#materials-science-community-forum-at-matsci.org)    Materials Science Community Forum at matsci.org

The Materials Project runs a forum at [matsci.org](https://matsci.org) intended as a shared space for several computational materials science projects, as well as general discussion about materials science. For the past several years, this effort has been co-run by the OpenKIM project. See [About matsci.org](https://matsci.org/faq) for more information about the forum and its governance.

All questions are welcome here! See our category at [https://matsci.org/materials-project](https://matsci.org/materials-project).

### [Direct link to heading](\#contact)    Contact

Please reach out to us on the [forum](https://matsci.org/materials-project): if questions or feedback are asked in a public setting, it allows others to benefit from seeing the answer too, and allows more people to participate in the conversation.

[PreviousDocumentation Credit](/credit) [NextGetting Involved](/community/getting-involved)

Last updated 2 years ago

---

# https://docs.materialsproject.org/community/getting-involved

The Materials Project would not be the resource it is today without the sustained efforts of many individual contributors who have helped make the Materials Project better. The Materials Project is a free, academic resource, with only a small team of core maintainers: any help received is always appreciated, and means we can make the Materials Project better for everyone!

There are several ways to get involved:

- If you are a software developer, please refer to the [Contributor Guide](/community/getting-involved/contributor-guide).

- If you are a domain expert, you can join the discussion and help answer questions of less experienced users in our forum at [https://matsci.org/materials-project](https://matsci.org/materials-project).

- If you are a domain expert, you can also notify us of errors, either in our public forum or via email at [feedback@materialsproject.org](mailto:feedback@materialsproject.org). Please check our [FAQ](/frequently-asked-questions) first to ensure that this error is not already known; some common issues arise from a misunderstanding of the data that Materials Project offers.

- If you generate data, either experimental or computational, you can use our contribution platform [MPContribs](/services/mpcontribs) to upload and link your data to the relevant material on Materials Project. This helps us by being able to offer a more complete and helpful resource, and also helps improve the discoverability of your own research by making it available to a wider audience. All uploaded data is credited to the original authors, and will have links to the appropriate publications.

- If you are an advanced user of Materials Project data or codes, you can help us improve documentation and tutorials.

- If you have discovered or know about a new crystal structure that is not present in the Materials Project database, you can submit it to us for calculation to help us offer a more complete database. If you are an advanced user, we may be able to receive calculations directly, but this typically requires prior communication and planning.


Any help is gratefully received, and we work hard to try to give back to the community ourselves wherever possible!

[PreviousGetting Help](/community/getting-help) [NextContributor Guide](/community/getting-involved/contributor-guide)

Last updated 11 months ago

---

# https://docs.materialsproject.org/community/getting-involved/contributor-guide

## [Direct link to heading](\#purpose-of-this-guide)    Purpose of this Guide

This guide aims to facilitate the process of contributing to any [Materials Project (MP) open source repositories](https://github.com/materialsproject). It offers high-level instructions and guidelines for those who wish to contribute to MP, regardless of the size of the contribution. Whether you're fixing a bug, improving docs, or proposing a new feature, this guide is for you. All contributions are welcome and appreciated!

This guide is a work in progress and will be updated as necessary to reflect changes in MP's practices. If you have any suggestions, don't hesitate to open an issue or a pull request!

**Happy contributing!**

## [Direct link to heading](\#the-materials-project-software-ecosystem)    The Materials Project Software Ecosystem

MP consists of several interconnected parts, each serving a specific purpose that together enable high-throughput computations.

### [Direct link to heading](\#official-materials-project-codes)    Official Materials Project Codes

The primary codes that most users are likely to interact with and contribute to are:

- **pymatgen**\[ [docs](https://pymatgen.org/)\]\[ [repo](https://github.com/materialsproject/pymatgen)\]: A large Python library for various materials analysis, manipulation and IO between different codes. Can be used on its own for analysis and setting up calculations to be executed manually, or together with the other codes below for a higher level of automation, error correction, and databasing of results.

- **atomate2**\[ [docs](https://materialsproject.github.io/atomate2/)\]\[ [repo](https://github.com/materialsproject/atomate2)\]: A library of automated computational workflows for various properties, such as structural relaxations, bandgaps, etc.


The following lower-level codes provide additional critical functions, but most users will likely make contributions to `pymatgen` or `atomate2`.

- **custodian**\[ [docs](https://materialsproject.github.io/custodian/)\]\[ [repo](https://github.com/materialsproject/custodian)\] **:** Just in time (JIT) job management software that provides automated, on the fly error correction to calculations as they are running. Many workflows (particularly VASP and Q-Chem) are designed to work closely with `custodian`, although use of custodian is not required.

- **jobflow** \[ [docs](https://materialsproject.github.io/jobflow/)\]\[ [repo](https://github.com/materialsproject/jobflow)\] **:** A library for writing and executing workflows. Jobflow defines the base `Job`, `Flow`, and `Maker` classes that are used in `atomate2` to define computational workflow steps.

- **fireworks**\[ [docs](https://materialsproject.github.io/fireworks/)\]\[ [repo](https://github.com/materialsproject/fireworks)\]: A software for managing execution of computational workflows, particularly suited for high-performance computing (HPC) environments with queueing systems. Instructions for setting up FireWorks for use with `atomate2` can be found [here](https://materialsproject.github.io/jobflow/install_fireworks.html). `atomate2` workflows can also be run without FireWorks.

- **emmet** \[ [docs](https://materialsproject.github.io/emmet/)\]\[ [repo](https://github.com/materialsproject/emmet/)\] **:** Defines structured schemas for storing outputs of different types of calculations performed by the Materials Project team. These comprise both code-specific schemas (e.g., for a VASP relaxation) and code-agnostic schemas (e.g., for any periodic solid material). `emmet` also uses maggma's `Builder` to define data processing pipelines that build the Materials Project database.

- **maggma**\[ [docs](https://materialsproject.github.io/maggma/)\]\[ [repo](https://github.com/materialsproject/maggma)\] **:** A framework for building modular data pipelines. maggma's `Store` and `Builder` classes provide a unified interface for accessing and transforming data. `atomate2` uses `Store` to save workflow results into a database or file, and `emmet` uses `Builder` to define the pipelines for processing Materials Project data.

- **crystaltoolkit** \[ [docs](https://docs.crystaltoolkit.org/index.html)\]\[ [repo](https://github.com/materialsproject/crystaltoolkit)\] **:** A web app framework that makes it easy for developers to create interactive web apps for materials science data, based on [plot.ly dash](https://dash.plotly.com/).


Because official MP codes are highly interdependent, their development is coordinated by the [MP Software Foundation](https://github.com/materialsproject/foundation). This group of developers meets regularly to establish policies regarding the scope of different packages, coding standards, etc.

### [Direct link to heading](\#external-tools-and-third-party-codes)    External Tools and Third-Party Codes

Many external or third-party developed codes are built to interoperate with the official MP codes above. An overview of these is available on the [Software Ecosystem page](/community/getting-involved/mp-community-software-ecosystem).

## [Direct link to heading](\#how-to-contribute)    How to Contribute

This section provides general guidelines for how to make contributions to the MP software ecosystem. If you are brand new to contributing to a software project, we encourage you to first read [Questions and Answers for new Contributors](/community/getting-involved/contributor-guide#questions-and-answers-for-new-contributors).

Note that detailed instructions for setting up a development environment or installing the necessary packages and dependencies for a particular project are not found here. Because they are repository-specific, please consult the documentation of the respective repositories (linked above) for that.

### [Direct link to heading](\#types-of-contributions)    Types of Contributions

We welcome many types of contributions, some of which require little to no coding experience. Contributions may include:

- Reporting a problem via a GitHub issue

- Testing a new feature

- Proposing a new feature

- Writing documenation

- Writing examples

- Developing graphics, slides, or Jupyter notebooks that aid in training and documentation

- Fixing a bug and submitting a GitHub pull request

- Writing a new feature and submitting a GitHub pull request


### [Direct link to heading](\#communication)    Communication

As you work on a contribution, the best ways to communicate with project maintainers and fellow users are:

- **GitHub Issues**: If you've found a problem or want to propose an idea, open an issue in the relevant repo. This is the first place to go if you need help with something. **Please don't submit how-to and support questions via issues, use GitHub Discussions instead (see below).**

- **Pull Request Comments**: If you want to discuss a specific change proposed in a pull request, use the PR's comments. This allows all discussions about a change to be kept in one place which is easily referenced later.

- **GitHub Discussions**: For more general discussions, use GitHub Discussions. This can be a great place to announce your intent to develop a new feature, ask for feedback on a proposal, discuss a new out-there idea, or get help with a problem.


Remember, it's okay to ask for help and feedback! We all started somewhere, and the MP community is there to help.

### [Direct link to heading](\#code-of-conduct)    Code of Conduct

Official Materials Project codes implement the [Contributor Covenant code of conduct](https://www.contributor-covenant.org/), which applies to project maintainers as well as all interactions with and among contributors. The overarching principle is to maintain a respectful and inclusive environment. Please read and adhere to these guidelines to ensure a positive and welcoming atmosphere for all contributors.

**TODO - need a link over "these guidelines" (need new PR to implement this)**

### [Direct link to heading](\#contribution-workflow)    \\#\# Contribution Workflow

Materials Project codes are hosted on GitHub and generally follow the [GitHub Flow](https://docs.github.com/en/get-started/quickstart/github-flow) development model. If you're unfamiliar with this process, refer to GitHub [docs](https://docs.github.com/en/github/collaborating-with-issues-and-pull-requests) for more information. Briefly, the steps are:

1. **Read this guide**: It provides an important orientation to the overall MP software ecosystem and expectations of code quality, etc.

2. **Check the Discussion boards:** Visit the GitHub discussion board for the project you want to contribute to to see whether anyone is working on something similar. You might find some free help!

3. **Describe your plans:** It's a good idea to post something on the discussion board to register your intentions (especially if you are developing a significant new feature or tutorial). This helps prevent duplication of effort.

4. **Fork and Clone**: Fork the repository you wish to contribute to, then clone it to your local machine.

5. **Create a New Branch**: Always create a new branch for your changes. This keeps your fork's main branch clean and makes it easier to open new pull requests in the future.

6. **Commit Your Changes**: Make your changes and commit them to your local repository.

7. **Push to Your Fork**: Push your changes to your fork.

8. **Open a Pull Request (PR)**: Open a PR against the **upstream** repo you're contributing to. We encourage you to do so EARLY - well before your code is highly developed or even working. You can use the Draft status to show that your PR is not ready for review yet, but having it open allows you to receive feedback from project maintainers and other community stakeholders. You can always mark it as ready for review later.


### [Direct link to heading](\#code-quality-guidelines)    \\#\# Code Quality Guidelines

Each Materials Project repository adheres to similar code format, testing, and documentation requirements. Although the specifics may vary slightly from repository to repository, the general requirements are as follows:

1. **Code Style**: All code should adhere to [`black`](https://black.readthedocs.io/en/stable/index.html) formatting and [`ruff`](https://beta.ruff.rs/docs/) linting rules, use `snake_case` for variable naming, `PascalCase` for classes, `CONSTANT_CASE` for globals.

2. **Testing**: All new features and bug fixes need tests. These should be implemented using [`pytest`](https://docs.pytest.org/en/7.3.x/), with a new unit test for each bug fix (that fails without the fix and passes with it) and functional tests for each new feature.

3. **Documentation**: Good docs are crucial. Function docstrings should follow the [Google docstring format](https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html) and describe every argument and keyword argument in concise terms, including appropriate units for the input (where applicable). Package documentation should be written in active and concise language, with small, ready-to-run code snippets that allow users to quickly try out new features. Relevant links should be included to allow users to easily find additional context or details e.g. in docs or GitHub issues/pull requests.


## [Direct link to heading](\#questions-and-answers-for-new-contributors)    Questions and Answers for new Contributors

### [Direct link to heading](\#how-much-experience-do-i-need-to-contribute)    How much experience do I need to contribute?

None! We routinely have new contributors who have not previously been involved in software development, or who are currently learning it as part of their graduate training. We do not have the resources to provide individual mentorship, but we will do what we can to support new contributors.

### [Direct link to heading](\#how-do-i-get-started)    How do I get started?

Reading this guide is the first step! After that, we suggest you visit the Discussion board of the GitHub repository for the project you want to contribute to. You can see what people are working on and even make a post to describe what you'd like to contribute to gather feedback.

If you prefer to discuss your plans more privately, don't be shy about reaching out to the package maintainers or other expert users directly.

### [Direct link to heading](\#where-can-i-find-help)    Where can I find help?

GitHub Issues for bugs and Discussions for Q&A (see [Communication](/community/getting-involved/contributor-guide#communication)) are a great place to start. For scientific and troubleshooting questions, you can also post on the [MatSci forums](https://matsci.org/). Finally, reach out to your colleagues, other expert users, or project maintainers.

### [Direct link to heading](\#how-do-i-add-a-new-workflow-link-to-atomate)    How do I add a new workflow? (link to atomate)

See this tutorial in the `atomate2` documentation! Note that all new workflows should go into `atomate2` rather than the legacy version of `atomate` (a.k.a., atomate 1).

**TODO - add link to atomate2 workflow tutorial**

### [Direct link to heading](\#how-do-i-support-a-new-code-in-pymatgen)    How do I support a new code in pymatgen?

See this tutorial in the `pymatgen` documentation. You can also draw inspiration from similar PRs. The most recent new code support was parsing AIRSS (ab-initio random structure search) results implemented in [pymatgen#2625](https://github.com/materialsproject/pymatgen/pull/2625).

**TODO develop tutorial and add link to pymatgen**

### [Direct link to heading](\#how-do-i-make-a-web-app-to-share-my-data)    How do I make a web app to share my data?

We suggest using `crystaltoolkit` , which we built to make it easy to create web apps for materials science. We have a [growing list of example apps on GitHub](https://github.com/materialsproject/crystaltoolkit/tree/main/crystal_toolkit/apps/examples) like this simple starter for rendering an [interactive 3d crystal structure](https://github.com/materialsproject/crystaltoolkit/blob/main/crystal_toolkit/apps/examples/basic_hello_structure_interactive.py). You can find a simple tutorial here.

**TODO - link to crystaltoolkit tutorial**

### [Direct link to heading](\#how-do-i-know-who-maintains-xxx-code)    How do I know who maintains XXX code?

Check the README file, which is displayed on the main page of each GitHub repository. We do our best to list the currently active maintainers of each repository there

**TODO - can/should we make this a policy?**

### [Direct link to heading](\#getting-credit-for-your-work)    Getting Credit for your Work

We value community contributions and want to do our best to provide appropriate credit. Some of the ways you can get visible credit for your work include

- Submit a PR to be added to the lists of contributors for a specific code. For example, see

- All MP codes support [`duecredit`](https://github.com/duecredit/duecredit/), which provides function decorators to associate publications with specific functions, classes and modules. You are welcome to include these decorators in your contributions where applicable (e.g. when you re-implement code from a paper, or use parameters from a paper, or contribute code you've created and published about). If in doubt, better to add citations than to not have them.

- If you have developed an external tool that uses one or more MP codes, we invite you to submit it for inclusion on the **TODO - link to** [Andrew Rosen](mailto:asrosen@lbl.gov) ecosystem page

- For `pymatgen` add-ons, submit a PR to be added to the [addons page](https://matsci.org/)


[PreviousGetting Involved](/community/getting-involved) [NextPotential Collaborators](/community/getting-involved/potential-collaborators)

Last updated 9 months ago

---

# https://docs.materialsproject.org/community/getting-involved/potential-collaborators



---

# https://docs.materialsproject.org/community/getting-involved/mp-community-software-ecosystem



---

# https://docs.materialsproject.org/community/community-resources



---

# https://docs.materialsproject.org/community/code-of-conduct



---

# https://docs.materialsproject.org/services/mpcontribs



---

# https://docs.materialsproject.org/methodology/materials-methodology



---

# https://docs.materialsproject.org/methodology/materials-methodology/overview



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/parameters-and-convergence



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/hubbard-u-values



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/gga+u-calculations/pseudopotentials



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/r2scan-calculations



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/r2scan-calculations/parameters-and-convergence



---

# https://docs.materialsproject.org/methodology/materials-methodology/calculation-details/r2scan-calculations/pseudopotentials



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/anion-and-gga-gga+u-mixing



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/thermodynamic-stability/gga-gga+u-r2scan-mixing



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/phase-diagrams-pds



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/chemical-potential-diagrams-cpds



---

# https://docs.materialsproject.org/methodology/materials-methodology/thermodynamic-stability/finite-temperature-estimation



---

# https://docs.materialsproject.org/methodology/materials-methodology/electronic-structure



---

# https://docs.materialsproject.org/methodology/materials-methodology/phonon-dispersion



---

# https://docs.materialsproject.org/methodology/materials-methodology/diffraction-patterns



---

# https://docs.materialsproject.org/methodology/materials-methodology/aqueous-stability-pourbaix



---

# https://docs.materialsproject.org/methodology/materials-methodology/magnetic-properties



---

# https://docs.materialsproject.org/methodology/materials-methodology/elasticity



---

# https://docs.materialsproject.org/methodology/materials-methodology/piezoelectric-constants



---

# https://docs.materialsproject.org/methodology/materials-methodology/dielectricity



---

# https://docs.materialsproject.org/methodology/materials-methodology/equations-of-state



---

# https://docs.materialsproject.org/methodology/materials-methodology/x-ray-absorption-spectra



---

# https://docs.materialsproject.org/methodology/materials-methodology/surface-energies



---

# https://docs.materialsproject.org/methodology/materials-methodology/grain-boundaries



---

# https://docs.materialsproject.org/methodology/materials-methodology/charge-density



---

# https://docs.materialsproject.org/methodology/materials-methodology/suggested-substrates



---

# https://docs.materialsproject.org/methodology/materials-methodology/related-materials



---

# https://docs.materialsproject.org/methodology/materials-methodology/optical-absorption-spectra



---

# https://docs.materialsproject.org/methodology/materials-methodology/alloys



---

# https://docs.materialsproject.org/methodology/molecules-methodology



---

# https://docs.materialsproject.org/methodology/molecules-methodology/overview



---

# https://docs.materialsproject.org/methodology/molecules-methodology/calculation-details



---

# https://docs.materialsproject.org/methodology/molecules-methodology/atomic-partial-charges



---

# https://docs.materialsproject.org/methodology/molecules-methodology/atomic-partial-spins



---

# https://docs.materialsproject.org/methodology/molecules-methodology/bonding



---

# https://docs.materialsproject.org/methodology/molecules-methodology/metal-coordination-and-binding



---

# https://docs.materialsproject.org/methodology/molecules-methodology/natural-atomic-and-molecular-orbitals



---

# https://docs.materialsproject.org/methodology/molecules-methodology/redox-and-electrochemical-properties



---

# https://docs.materialsproject.org/methodology/molecules-methodology/molecular-thermodynamics



---

# https://docs.materialsproject.org/methodology/molecules-methodology/vibrational-properties



---

# https://docs.materialsproject.org/methodology/molecules-methodology/legacy-data



---

# https://docs.materialsproject.org/methodology/mof-methodology



---

# https://docs.materialsproject.org/methodology/mof-methodology/calculation-parameters



---

# https://docs.materialsproject.org/methodology/mof-methodology/calculation-parameters/dft-parameters



---

# https://docs.materialsproject.org/methodology/mof-methodology/calculation-parameters/density-functionals



---

# https://docs.materialsproject.org/methodology/mof-methodology/calculation-parameters/psuedopotentials



---

# https://docs.materialsproject.org/methodology/mof-methodology/calculation-parameters/dft-workflow



---

# https://docs.materialsproject.org/apps/explorer-apps



---

# https://docs.materialsproject.org/apps/explorer-apps/materials-explorer



---

# https://docs.materialsproject.org/apps/explorer-apps/materials-explorer/tutorial



---

# https://docs.materialsproject.org/apps/explorer-apps/molecules-explorer



---

# https://docs.materialsproject.org/apps/explorer-apps/molecules-explorer/tutorial



---

# https://docs.materialsproject.org/apps/explorer-apps/molecules-explorer/legacy-data



---

# https://docs.materialsproject.org/apps/explorer-apps/battery-explorer



---

# https://docs.materialsproject.org/apps/explorer-apps/battery-explorer/background



---

# https://docs.materialsproject.org/apps/explorer-apps/battery-explorer/tutorial



---

# https://docs.materialsproject.org/apps/explorer-apps/synthesis-explorer



---

# https://docs.materialsproject.org/apps/explorer-apps/synthesis-explorer/background



---

# https://docs.materialsproject.org/apps/explorer-apps/synthesis-explorer/tutorial



---

# https://docs.materialsproject.org/apps/explorer-apps/catalysis-explorer



---

# https://docs.materialsproject.org/apps/explorer-apps/catalysis-explorer/tutorial



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/downloading-the-data



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/structure-details



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/structure-details/qmof-ids



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/structure-details/structure-sources



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/structure-details/finding-mofs-by-common-name



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/structure-details/structural-fidelity



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions/smiles-mofid-and-mofkey



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions/pore-geometry



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions/topology



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions/electronic-structure



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions/population-analyses-and-bond-orders



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/property-definitions/symmetry



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/version-history



---

# https://docs.materialsproject.org/apps/explorer-apps/mof-explorer/how-to-cite



---

# https://docs.materialsproject.org/apps/analysis-apps



---

# https://docs.materialsproject.org/apps/analysis-apps/phase-diagram



---

# https://docs.materialsproject.org/apps/analysis-apps/phase-diagram/background



---

# https://docs.materialsproject.org/apps/analysis-apps/phase-diagram/tutorials



---

# https://docs.materialsproject.org/apps/analysis-apps/phase-diagram/faq



---

# https://docs.materialsproject.org/apps/analysis-apps/pourbaix-app



---

# https://docs.materialsproject.org/apps/analysis-apps/pourbaix-app/background



---

# https://docs.materialsproject.org/apps/analysis-apps/pourbaix-app/tutorial



---

# https://docs.materialsproject.org/apps/analysis-apps/pourbaix-app/faq-+-known-issues



---

# https://docs.materialsproject.org/apps/analysis-apps/crystal-toolkit



---

# https://docs.materialsproject.org/apps/analysis-apps/crystal-toolkit/background



---

# https://docs.materialsproject.org/apps/analysis-apps/crystal-toolkit/tutorial



---

# https://docs.materialsproject.org/apps/analysis-apps/crystal-toolkit/faq



---

# https://docs.materialsproject.org/apps/analysis-apps/reaction-calculator



---

# https://docs.materialsproject.org/apps/analysis-apps/interface-reactions



---

# https://docs.materialsproject.org/apps/characterization-apps



---

# https://docs.materialsproject.org/apps/characterization-apps/x-ray-absorption-spectra-xas



---

# https://docs.materialsproject.org/apps/explore-contributed-data



---

# https://docs.materialsproject.org/downloading-data/how-do-i-download-the-materials-project-database



---

# https://docs.materialsproject.org/downloading-data/using-the-api



---

# https://docs.materialsproject.org/downloading-data/using-the-api/getting-started



---

# https://docs.materialsproject.org/downloading-data/using-the-api/querying-data



---

# https://docs.materialsproject.org/downloading-data/using-the-api/tips-for-large-downloads



---

# https://docs.materialsproject.org/downloading-data/using-the-api/examples



---

# https://docs.materialsproject.org/downloading-data/using-the-api/advanced-usage



---

# https://docs.materialsproject.org/downloading-data/using-the-api/legacy-api-users



---

# https://docs.materialsproject.org/downloading-data/using-the-api/deprecated-materials

The Materials Project is constantly updating its data with more accurate calculations and improved standards. As a result, some materials are "deprecated" to facilitate access to the highest quality data available, which is where older calculations are hidden from default search results but are still searchable using our API and website. A list of deprecated entries and their reasons for deprecation is maintained in the [Materials Project Github user guide](https://github.com/materialsproject/docs/blob/master/docs/user-guide/deprecations.md).

[PreviousLegacy API Users](/downloading-data/using-the-api/legacy-api-users) [NextDifferences between new and legacy API](/downloading-data/differences-between-new-and-legacy-api)

Last updated 2 years ago

---

# https://docs.materialsproject.org/downloading-data/differences-between-new-and-legacy-api

A summary of differences between the new and legacy API can be found in the table below. For more detailed information on differences regarding the Python API client, please see the [Legacy API Users](/downloading-data/using-the-api/legacy-api-users) section.

|  | New API | Legacy API |
| --- | --- | --- |
| **Currently recommended for** | Early adopters | Everyone else |
| **Base URL** | api.materialsproject.org | materialsproject.org/rest/v2 |
| **Documentation** | [api.materialsproject.org/docs](https://api.materialsproject.org/docs) | [mapidoc](https://github.com/materialsproject/mapidoc) |
| **Specification** | [OpenAPI-compliant specification available](https://api.materialsproject.org/openapi.json) | None available |
| **Support** | Our new API will be supported for the forseeable future once released | Will be available for at least one year after new API is finalized |
| **Data Updates** | Will receive new data updates included latest and most accurate data | Will be frozen at database release v2021.03.13 |
| **API Key** | [Available here](https://next-gen.materialsproject.org/api#api-key) | Available at [legacy.materialsproject.org/open](https://legacy.materialsproject.org/open) |
| **Python client installation** | `pip install mp-api` (may be available in _pymatgen_ at a later date) | `pip install pymatgen` |
| **Python client import code** | `from mp_api.client import MPRester` | `from pymatgen.ext.matproj import MPRester` |
| **MPContribs integration for user contributed data** | Yes | No |

[PreviousDeprecated Materials](/downloading-data/using-the-api/deprecated-materials) [NextQuery and Download Contributed Data](/downloading-data/query-and-download-contributed-data)

Last updated 2 years ago

---

# https://docs.materialsproject.org/downloading-data/query-and-download-contributed-data

## [Direct link to heading](\#interactively)    Interactively

Under Construction. We are working on enabling CSV/JSON download buttons for entire projects or specific queries. The downloads will also function as versioned snapshots for each project.

## [Direct link to heading](\#programmatically)    Programmatically

Querying MPContribs data programmatically involves the following steps:

- Install the python package [mpcontribs-client](https://pypi.org/project/mpcontribs-client/).

- Initialize the client with your API key and a project

- Check `client.available_query_params()` for available query parameters

- Build up a query dictionary

- Decide which fields to retrieve

- Run `client.query_contributions()`


Here's an example for the `carrier_transport` dataset:

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
from mpcontribs.client import Client
client = Client(apikey="your-api-key-here", project="carrier_transport")
client.available_query_params()  # print list of available query parameters
query = {"formula__contains": "Au", "data__PF__p__value__lt": 10}
fields = ["identifier", "formula", "data.metal", "data.S.n.value"]
data = client.query_contributions(
    query=query, fields=fields, sort="-data.S.n.value", paginate=True
)
```

By default, `paginate` is `False` which will only retrieve the first page of results and should be used to test the `query`, `fields` and `sort` parameters before paginating through all results.

If entire projects or large subsets of contributed data are downloaded for later use, it is often more efficient to use the `client.download_contributions()` function. It also takes a `query` as argument and downloads all results as `json.gz` files behind the scenes. Only locally missing data is downloaded when `download_contributions` is run and contributions are loaded from disk. This function always retrieves all fields included in the `data` component, so the `fields` argument is not available/needed. Additional components (i.e. `structures`, `tables`, and `attachments`) can be included in the downloads through the `include` argument:

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
data = client.download_contributions(query=query, include=["tables"])
```

[PreviousDifferences between new and legacy API](/downloading-data/differences-between-new-and-legacy-api) [NextTroubleshooting API Issues](/downloading-data/troubleshooting-api-issues)

Last updated 7 months ago

---

# https://docs.materialsproject.org/downloading-data/troubleshooting-api-issues

If you have trouble using the Materials Project API, please try the following before asking questions on our forum. If you're still having problems, make sure to report your Python version, `mp-api` package version, the specific line of code you're trying to run, and any other relevant information about your system.

- Try updating the `mp-api` package by `pip install --upgrade mp-api` or similar.


[PreviousQuery and Download Contributed Data](/downloading-data/query-and-download-contributed-data) [NextAWS OpenData](/downloading-data/aws-opendata)

Last updated 1 year ago

---

# https://docs.materialsproject.org/downloading-data/aws-opendata

In an effort to make our data as accessible as possible (FAIR principle) as well as significantly improve data downloads and take pressure off our servers, we are making a growing list of our data products available through the [AWS OpenData Program](https://aws.amazon.com/opendata). Also see the entries for MP-managed data on the [AWS OpenData Registry](https://registry.opendata.aws/materials-project/) or the [AWS Data Exchange](https://aws.amazon.com/marketplace/pp/prodview-hc3sh3u5ukiya).

## [Direct link to heading](\#overview)    Overview

MP data is organized in 3 buckets named `materialsproject-{raw,parsed,build}`. Note that the particular organization of our data in these buckets is still in flux and can change without notice as we integrate them into our cloud infrastructure.

### [Direct link to heading](\#raw-data)    raw data

We are in the process of providing VASP output files for our calculations in the `raw` bucket. Look out for announcements through our email lists and notifications on our website.

### [Direct link to heading](\#parsed-data)    parsed data

The [`parsed` bucket](https://materialsproject-parsed.s3.amazonaws.com/index.html) contains objects that MP generates by parsing the VASP output files. The objects form the basis for our builder pipelines which create the derived high-level data collections served through the MP API and website. All S3 objects in this bucket are serialized `pymatgen` or `emmet` python objects and most are stored as gzip-compressed JSON files for each MP ID (i.e. `<prefix>/<mp-id>.json.gz`). We are in the process of grouping documents into JSON Lines (JSONL) files to reduce the number of files and significantly improve transfer speeds. `tasks` are now organized by `nelements/output.spacegroup.number` and a timestamp ( `dt`) derived from the earliest `completed_at` in the list of tasks included in the respective object.

| prefix | \# objects | size |
| --- | --- | --- |
| `/dos` | 691k | 63.1 GB |
| `/bandstructures` | 705k | 1.4 TB |
| `/chgcars` | 400k | 7.2 TB |
| `/aeccar{0,1,2}s` | 138.7k each | 1.1 TB each |
| `/elfcars` | 107.5k | 101 GB |
| `/locpots` | 158k | 2.5 TB |
| `/tasks` | 1556 | 34 GB |

### [Direct link to heading](\#build-data)    build data

The [`build` bucket](https://materialsproject-build.s3.amazonaws.com/index.html) contains the high-level derived data that comprises the source for the collections available through the [MP API](https://api.materialsproject.org/) as well as pre-built objects and images for efficient visualization on the website.

The collections and pre-built objects are versioned by the database release date and individual documents grouped into gzip-compressed JSONL files. Images are stored in PNG format. Use the `ls` command for the AWS CLI or the [bucket explorer](https://materialsproject-build.s3.amazonaws.com/index.html) to list the categories available under each prefix (see [download section](/downloading-data/aws-opendata#explore-and-download-data) below).

| prefix | version | \# objects | size |
| --- | --- | --- | --- |
| `/collections` | `/2022-10-28` | 12.6k | 2.8 GB |
|  | `/2023-11-01` | 18.4k | 6.1 GB |
| `/objects` | `/2022-10-28` | 289k | 55.9 GB |
| `/images` | N/A | 200k | 58 GB |

## [Direct link to heading](\#explore-and-download-data)    Explore & Download Data

We are in the process of integrating all available data into the `mp-api` python client for improved convenience and efficiency. However, all data in MP's OpenData buckets can always be downloaded directly using the [AWS CLI](https://docs.aws.amazon.com/cli/latest/userguide/getting-started-install.html).

Start by exploring the contents of the bucket you're interested in, by either navigating to the bucket's web interface (e.g. [https://materialsproject-parsed.s3.amazonaws.com/index.html](https://materialsproject-parsed.s3.amazonaws.com/index.html)) or using the CLI's `ls` command:

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
aws s3 ls --no-sign-request s3://<bucket>/<prefix>/[<version>]/

# examples
aws s3 ls --no-sign-request s3://materialsproject-parsed/
aws s3 ls --no-sign-request s3://materialsproject-build/
aws s3 ls --no-sign-request s3://materialsproject-build/collections/2022-10-28/
aws s3 ls --no-sign-request s3://materialsproject-build/images/
```

All objects for a prefix can be downloaded, using the format

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
# the AWS CLI will parallelize download requests automatically
aws s3 cp --no-sign-request --recursive s3://<bucket>/<prefix>/ <output-dir>/

# examples
aws s3 cp --no-sign-request --recursive s3://materialsproject-parsed/tasks/ mp-tasks/
aws s3 cp --no-sign-request --recursive \
    s3://materialsproject-build/images/structures/ mp-images-structures/
```

[PreviousTroubleshooting API Issues](/downloading-data/troubleshooting-api-issues) [NextContribute Data](/uploading-data/what-is-mpcontribs)

Last updated 7 months ago

---

# https://docs.materialsproject.org/uploading-data/what-is-mpcontribs

### [Direct link to heading](\#interactively)    Interactively

Under Construction. We are working on creating interactive interfaces that enable our user to contribute data in addition to the programmatic approach below.

### [Direct link to heading](\#programmatically)    Programmatically

#### [Direct link to heading](\#quickstart)    Quickstart

Read the [concepts section](/services/mpcontribs#concepts), [Create a project](https://contribs.materialsproject.org/contribute), install [`mpcontribs-client`](https://pypi.org/project/mpcontribs-client/), and set the `MPCONTRIBS_API_KEY` environment variable to the API key shown on the [profile page](https://profile.materialsproject.org) or your [dashboard](https://materialsproject.org/dashboard). The following code snippet outlines the general process of adding data to your project. See the next section for step-by-step instructions.

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
from mpcontribs.client import Client
client = Client(project="my_test_project")
columns = {"a": "eV", "b.c": None, "b.d": ""}
client.init_columns(columns)
contributions = [{\
    "identifier": "mp-4",\
    "data": {\
        "a": "3 eV",\
        "b": {"c": "hello", "d": 3}\
    },\
    "structures": [`pymatgen.core.Structure`, ...],\
    "tables": [`pandas.DataFrame`, ...],\
    "attachments": [`pathlib.Path`, `mpcontribs.client.Attachment`, ...]\
}]
client.submit_contributions(contributions)
client.init_columns(columns) # shouldn't be needed but ensures all columns appear
```

#### [Direct link to heading](\#step-by-step-instructions)    Step-by-step instructions

The general process of preparing and contributing data using the [mpcontribs-client python library](https://pypi.org/project/mpcontribs-client/) follows the steps below. Make sure to read the [Concepts section](/services/mpcontribs#concepts) before continuing.

**Install dependencies.** You might not need all of them depending on your use case:
`pip install mpcontribs-client mp_api pandas flatten_dict pymatgen monty`

**Import commonly used libraries** (some might be optional in your case):

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
import json
import os
import gzip

import pandas as pd
import numpy as np

from glob import glob
from pathlib import Path
from math import isnan
from flatten_dict import flatten, unflatten
from collections import defaultdict

from monty.serialization import loadfn
from pymatgen.core import Structure
from mp_api.client import MPRester
from mpcontribs.client import Client, Attachments
```

[Initialize `MPRester`](/downloading-data/using-the-api/getting-started) and **create a project.** Take an extra minute to think about the `name` for your project. It needs to be 3-31 characters long and can only use alpha-numeric characters or the underscore ( `_`). You'll use it to refer to your project in the python client, and it'll be used as part of the URL for your project's landing page.

You can either create the project by [completing this form](https://contribs.materialsproject.org/contribute) or programmatically as shown below. Replace the example values below with the info pertinent to your project. Projects created with these example values or with insufficient information will be removed immediately. Also check out [the doc strings](https://jakevdp.github.io/PythonDataScienceHandbook/01.01-help-and-documentation.html) for more info about any functions used here.

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
with MPRester() as mpr:  # needs MP_API_KEY environment variable to be set
    mpr.contribs.create_project(
        name="my_test_project",
        title="My Test Title",
        authors="P. Huck, K. Persson",
        description="Describe the contents of your project here.",
        url="https://example.com/projects/test"
        # URL could be group website, publication DOI, Github repo, ...
    )
```

After the project is created, you can use the MPContribs client directly to exclusively interact with your project going forward:

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
# set MPCONTRIBS_API_KEY environment variable first (same as MP_API_KEY)
client = Client(project="my_test_project")
client.get_project()
```

**Update submitted project information** if needed. For instance, we might want to add another first author to `authors` or change the label for the default reference and add another URL to `references`:

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
client.update_project({
    "authors": "P. Huck / J. Munro, K. Persson",
    "references": [\
        {"label": "Website", "url": "https://example.com/projects/test"},\
        {"label": "PRL", "url": "https://doi.org/123456"}\
    ]
})
client.get_project(fields=["name", "authors", "references"])
```

**Load the data** you plan to contribute from disk. As an example,

- we load the following CSV file ( `main.csv`) with each line containing the values intended for the queryable columns in the `data` component for a material. Use `mpr.find_structure()` and/or `mpr.summary.search()` to identify MP materials matching your structure, if possible. Only contributions with MP IDs as `identifier` will show up on MP's materials details pages.







Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
identifier,energy,polarization,max_distance,polar_icsd,polar_spacegroup
mp-4,3.5,2,4,5.6,7
mp-100,15,33,22.4,16,8
```

- we load a list of associated spectra intended for the `tables` component from a folder name `spectra` with the following contents:







Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
spectra/
|- mp-4.csv
|    energy,XAS,XMCD
|    1,3,5
|    2,4,6
|- mp-100.csv
|    [...]
```

- we load a list of associated CIF files intended for the `structures` component as [pymatgen structures](https://pymatgen.org/usage.html) (folder named `structures` contains a list of CIF files named by MP ID).

- we load data intended for the `attachments` component directly from files/images on disk using the following directory structure. Attachments can also be created dynamically from standard python lists or dictionaries. See the use of `Attachments.from_list()` below or its doc string for more details.







Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
attachments/
|- mp-4/
|  |- image.jpg
|  |- extra.json
|  |- info.txt
|- mp-100/
|  [...]
```


Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
main = pd.read_csv("main.csv", index_col="identifier")
tables = {  # map mp-id to pandas DataFrame
    path.stem: pd.read_csv(path, index_col="energy")
    for path in Path("spectra").glob("*.csv")
}
structures = {
    path.stem: Structure.from_file(path)
    for path in Path("structures").glob("*.cif")
}
attachments = defaultdict(list)
for path in Path("attachments").rglob("*"):
    if path.is_file():
        attachments[path.parent.name].append(path)
```

**Initialize the columns** for the queryable `data` component of your contributions. This involves deciding on MPContribs-compatible field names, appropriately grouping / categorizing related columns, setting their units, and adding descriptions to be included in the project info.

- Use dot ( `.`) notation to group/nest columns up to 4 levels deep. Most non-alphanumeric characters are disallowed (including underscores and spaces) to encourage better data organization by grouping columns and using short, readable, and type-able column names. Use a good description to explain column names, and use the pipe ( `|`) character to indicate conditions for a column (e.g. `max`, `300K`, ...) where nesting might not be desired. Falling back on CamelCase is also an option for column naming.

- The column unit can either be an empty string ( `""`) to indicate a dimensionless number, a string representing a unit supported by [pint](https://pint.readthedocs.io/), or `None` to indicate that the column values are not numerical.

- To remove keys from the `other` field, explicitly set those keys to `None` in `update_project()`


Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
columns_map = {
    "energy": {"name": "energy", "unit":"eV/atom", "description": "Formation Energy"},
    "polarization": {"name": "polarization", "unit": "µC/cm²", "description": "..."}},
    "max_distance": {"name": "distance|max", "unit": "Å", "description": "..."},
    "polar_icsd": {"name": "polar.icsd", "unit": "", "description": "..."},
    "polar_spacegroup": {"name": "polar.spacegroup", "description": "..."}}
}
columns = {col["name"]: col.get("unit") for col in columns_map.values()}
client.init_columns(columns)
other = unflatten({
    col["name"]: col["description"] for col in columns_map.values()
}, splitter="dot")
client.update_project({"other": other})
```

**Prepare the list of contributions.**

- The preferred `data` format is to provide numerical data as strings with units rather than naked numerical values. This will be parsed by the API to yield a value and a unit for display and promotes the correct reporting of significant figures.

- If you're trying to include `list` objects in the `data` component, first convert them into dictionaries by giving each element in the list a name. For instance, in the case of tensors, convert `[[1, 2], [3, 4]]` to `{"e11": 1, "e12": 2, "e21": 3, "e22": 4}`. This will make the tensor components queryable.


Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
contributions = []

for record in main.to_dict("records"):
    clean = {}
    for k, v in record.items():
        # remove NaNs (tip: skip any unset/empty keys)
        if not isinstance(v, str) and isnan(v):
            continue
        # convert boolean values to Yes/No, and append units
        key = columns_map[k]["name"]
        unit = columns_map[k].get("unit")
        val = v
        if isinstance(v, bool):
            val = "Yes" if v else "No"
        elif unit:
            val = f"{v} {unit}"

        clean[key] = val

    identifier = clean.pop("identifier")
    contrib = {"identifier": identifier}
    contrib["data"] = unflatten(clean, splitter="dot")

    table = tables.get(identifier)
    if table:
        # set table attributes for plotting (see Plotly Express)
        table.attrs = {
            "name": "absorption coefficients",
            "title": "Energy-dependent Absorption Coefficients",
            "labels": {
                "value": "absorption coefficient [cm⁻¹]",
                "variable": "method"
            }
        }
        contrib["tables"] = [table]

    structure = structures.get(identifier)
    if structure:
        structure.name = "optional structure name"
        contrib["structures"] = [structure]

    attms = attachments.get(identifier)
    if attms:
        contrib["attachments"].from_list(attms)

    contributions.append(contrib)
```

**Submit the contributions and publish the project** when ready.

Copy

```min-w-full inline-grid [grid-template-columns:auto_1fr] py-2 px-2 [counter-reset:line]
client.submit_contributions(contributions)
client.init_columns(columns) # shouldn't be needed but ensures all columns appear
# client.delete_contributions()  # remove all contributions from project
client.make_public()
```

### [Direct link to heading](\#api-docs-page)    API Docs Page

The [mpcontribs-client python library](https://pypi.org/project/mpcontribs-client/) interacts with the [MPContribs API](https://contribs-api.materialsproject.org) to programmatically access or retrieve experimental and theoretical data contributed by the MP community. Project information is retrievable through the `projects` resource, and the corresponding contributed data through the `contributions` resource. Each project can contain many contributions for an MP material or composition. Each contribution in turn consists of four (optional) components: free-form hierarchical data, tabular data, crystal structures, and attachments. There are separate dedicated resource endpoints for `tables`, `structures`, and `attachments`. See the [Concepts section](/services/mpcontribs#concepts) for more details. Check out the "Models" section on the [API Docs](https://contribs-api.materialsproject.org) page for descriptions of available fields and to try out the API in the browser.

[PreviousAWS OpenData](/downloading-data/aws-opendata) [NextData Workflows](/data-production/data-workflows)

Last updated 27 days ago

---

# https://docs.materialsproject.org/data-production/data-workflows

Materials Project data is produced using automated high-throughput workflows implemented within the open-source [atomate](https://atomate.org/) software package. A list of the workflows used to produce the core set of materials data include:

- **Static (wf\_static)**

- **Structure Optimization (wf\_structure\_optimization)**

- **Band Structure and Density of States (wf\_bandstructure)**

- **Bulk Modulus and Equations of State (wf\_bulk\_modulus)**

- **Dielectric Constants (wf\_dielectric\_constant)**

- **Elastic Constants (wf\_elastic\_constant)**

- **Piezoelectric Constants (wf\_piezoelectric\_constant)**

- **r2SCAN Optimization (wf\_r2scan\_opt)**

- **SCAN Optimization (wf\_scan\_opt)**

- **X-ray Absorption Spectra (wf\_xas, wf\_exafs\_paths, wf\_eels)**

- **Surface Energies (wf\_slab)**


The DFT calculations run by these workflows use standardized input sets [defined within pymatgen](https://pymatgen.org/pymatgen.io.vasp.sets.html). The goal of these workflows is to produce calculation (task) data which is parsed from the raw output files using a "drone" (e.g. [atomate's vasp drone](https://atomate.org/atomate.vasp.html?highlight=drone#module-atomate.vasp.drones)) and stored within a MongoDB collection. From this, other database collections are "built" to house specific data (i.e. thermo, dielectric, xas, etc...). For more information on builders and building data, see the [Data Builders](/data-production/data-builders) section.

[PreviousContribute Data](/uploading-data/what-is-mpcontribs) [NextData Builders](/data-production/data-builders)

Last updated 2 years ago

---

# https://docs.materialsproject.org/data-production/data-builders

Builders to produce MongoDB collections associated with specific categories of materials data are implemented in the [emmet-builders](https://github.com/materialsproject/emmet/tree/main/emmet-builders) software package. These take a set of MongoDB collections as input, extract and transform data from them, and then output a new collection. The schema used by builders for the input and output collections is defined by a set of standardized document models. We use a Python library called [pydantic](https://pydantic-docs.helpmanual.io) to structure these document models, and we store our documents models within the [emmet-core](https://github.com/materialsproject/emmet/tree/main/emmet-core) software package. Additionally, these models are used to define the schema for the Materials Project API which has its server-side code implemented in [emmet-api](https://github.com/materialsproject/emmet/tree/main/emmet-api).

To browse our document models defined in emmet, see [here](https://github.com/materialsproject/emmet/tree/main/emmet-core/emmet/core). For example, see the `ThermoDoc` defined [here](https://github.com/materialsproject/emmet/blob/main/emmet-core/emmet/core/thermo.py) as an example document model that powers both the [`ThermoBuilder`](https://github.com/materialsproject/emmet/blob/main/emmet-builders/emmet/builders/vasp/thermo.py) and the Materials Project's [thermo API endpoint](https://api.materialsproject.org/docs#/Thermo/get_by_key_thermo__material_id___get).

The figure below illustrates the entire Materials Project build pipeline including builders and all input/output collections:

![](https://2369879881-files.gitbook.io/~/files/v0/b/gitbook-x-prod.appspot.com/o/spaces%2F-MhdHkeirg8PPHDHWitE%2Fuploads%2FFB5MjsBqjkqiRG4byUDM%2Fmp_build_pipeline.svg?alt=media&token=aadf10db-388a-4f95-99c8-4ca21d77d109)

For information on how to run any of the emmet builders, see the [Running Builders](https://materialsproject.github.io/maggma/getting_started/running_builders/) section of the maggma software package which defines a lot of the core builder related code and CLI.

[PreviousData Workflows](/data-production/data-workflows)

Last updated 2 years ago

---

