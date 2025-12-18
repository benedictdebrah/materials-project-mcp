---
abstract: |
  The Materials Project (MP) has transformed computational materials
  science by providing open-access data on crystal structures,
  electronic properties, and thermodynamic stability. However, its
  dependence on programmatic queries and complex JSON responses creates
  accessibility barriers for non-expert users. We present the Materials
  Model Context Protocol (MCP) server, an MCP server integrating the MP
  API with large language models (LLMs) to enable natural language
  access to structured materials data. The server employs schema-driven
  tool definitions with rigorous input/output validation, asynchronous
  processing, and compatibility with MCP-aware clients such as Claude
  Desktop and VS Code Copilot. Case studies on LiFeO and BaTiO systems
  show accurate, context-rich data retrieval comparable to direct API
  queries, offering interactions that are more intuitive and adaptable.
  Evaluation confirms the MCP server preserves scientific accuracy,
  supports complex searches, and reduces the technical burden. This
  study highlights MCP-based interfaces' potential to democratize
  materials data access, accelerate hypothesis generation, and enable
  data-driven discovery.\
  **Keywords:** Materials Project, Model Context Protocol, Material MCP,
  API, Server
author:
- Benedict O. Debrah
- Etor P. Fiawornu
- Luke E. K. Achenie
bibliography:
- Citations.bib
title:  The Materials Project Platform - an MCP based server
---

# Introduction

The Materials Project (MP) platform, a pioneering repository for
computational materials science data, has improved the design and
discovery of novel materials by providing access to data associated with
properties, inorganic compounds and structure images[@osti_1511347].
However, its vast and heterogeneous datasets spanning crystal
structures, electronic properties and synthesis protocols become very
difficult for non-expert users due to the platform's reliance on
programmatic queries (e.g., Python APIs) and fragmented documentation.
Current accessibility of materials data from the Materials Project
platform is limited, primarily due to the need for programming skills
and the complexity of interpreting nested JSON responses from manual API
queries [@Jain2018] [@persson_design_nodate]. Accessibility restrictions
create a barrier that makes the data from the materials less accessible
to a wider range of researchers. Specifically, researchers face
challenges when attempting to write and execute complex queries for
information or data retrieval[@Himanen2019DataDriven].

To mitigate such limitations in scientific data access, Our Materials
Project Model Context Protocol (MCP) offers a standardized approach for
integrating data sources like the Materials Project database with large
language models (LLMs). Introduced by Anthropic in late
2024,[@modelcontextprotocol2025] MCP employs a client-server
architecture to streamline connections between LLMs and domain-specific
tools. The protocol has demonstrated success in lowering technical
barriers for database systems (e.g., MySQL, PostgreSQL), enabling
natural language queries for complex data retrieval. This paradigm shows
significant promise for democratizing access to specialized scientific
repositories.

## Materials Project

The Materials Project is a open-access database for computational
materials that has revolutionized materials science by offering an
extensive collection of computed properties for numerous inorganic
compounds [@Jain2018]. Established in 2011, the platform's primary aim
is to enhance materials discovery and design via high-throughput Density
Functional Theory (DFT) computations[@JAIN20112295]. It provides
extensive information, encompassing crystal structures, electronic band
gaps, and formation energies, essential for researchers looking to
forecast and evaluate new materials with specific functionalities,
including those for energy storage and conversion
[@horton2025accelerated].The platform's accessibility, via its web
interface and API, supports a data-centric approach to materials
research, in line with the tenets of the Materials Genome Initiative.
The Materials Project has cultivated a cooperative community by offering
its data and software infrastructure as open-source, enabling ongoing
growth and the creation of new analytical tools [@D3DD00153A]. Its
effect is shown through its extensive application, as hundreds of
thousands of researchers utilize its data to comprehend, forecast, and
create innovative materials[@becker2022progress].

## Other Related Work

Artificial intelligence (AI) and machine learning (ML) are fundamentally
transforming materials science, enabling researchers to generate more
breakthrough ideas and predict material properties before physical
synthesis [@jaiswal2024]. Although early applications of natural
language processing (NLP) in material science focused on extracting
information or summarizing crystal structures[@Olivetti2020DataDriven],
and prompt-based methods automated data
extraction[@chiang2024llamplargelanguagemodel]. Through these advances
in AI tools, advanced techniques emerged, such as AI agents, which have
helped automate materials discovery
[@Pasupuleti2024AIPowered; @Ghafarollahi2025Automating]. Knowledge
graphs have also advanced materials science by representing information
in structured networks that facilitate data accessibility and material
discovery[@venugopal2022matkglargestknowledgegraph]. One such Knowledge
Graph is MatKG, a prominent Knowledge Graph and the largest in materials
science, which illustrates this advancement by organizing a vast array
of materials science data in a structured
format[@venugopal2022matkglargestknowledgegraph]. Materials Knowledge
Graph (MKG), which uses NLP and LLM to extract and research structured
triples, and Materials Design Ontology (MDO), which aims to standardize
access to the computational materials
database[@ye2025constructionapplicationmaterialsknowledge]. There is an
automated Retrieval-Augmented Generation (RAG) framework proposed for
materials knowledge retrieval, which employs reasoning and acting
(ReAct) agents to interact with computational and experimental data from
the Materials Project; LLaMP can process higher order data, streamline
complex tasks in computational materials.
[@chiang2024llamplargelanguagemodel].

Despite these significant advances, current approaches sometimes face
several limitations. Artificial intelligence agents, while powerful for
research and automation, are often susceptible to
hallucinations[@bekahhw_agents_vs_mcps_2025]. Traditional Knowledge
Graphs (KGs), while successful, are often static and time-consuming to
generate output from[@jiang2025evolution], which struggle to keep up
with the rapidly evolving materials knowledge landscape. In addition,
these KGs require human expertise for annotation and ontology
development[@osti_10393958], which limits their ability to incorporate
new discoveries. These limitations present a critical gap in current
techniques, especially for the accessibility of materials data.

To address these challenges, the Model Context Protocol (MCP) has
recently emerged as a promising solution. MCP supports dynamic data
retrieval, making it particularly suited for domains characterized by
rapidly evolving knowledge, such as materials science.

In Section [2](#sec:methodology){reference-type="ref"
reference="sec:methodology"}, we will discuss the methodologies employed
in the development and implementation of the Materials Project MCP. In
Section [3](#sec:evaluation){reference-type="ref"
reference="sec:evaluation"}, we will evaluate the effectiveness and
performance of the Materials Project MCP. This includes usability and
case studies. Finally, in Section [4](#sec:summary){reference-type="ref"
reference="sec:summary"}, we will summarize the key findings and suggest
directions for future improvements.

# Methodology {#sec:methodology}

## System Architecture

This study aims to significantly enhance materials data accessibility by
designing and implementing a system that integrates AI applications via
the interface of the Model Context Protocol (MCP) and the materials
project API. This would enable programmatic and schema-driven access to
materials data through large language models (LLMs) and agentic
interfaces such as Claude Desktop, VScode Copilot [@Narajala2025MCP].
Our Materials MCP server connects natural language requests from users
to the organized API endpoints offered by the Materials Project.

### Schematic Overview {#schematic-overview .unnumbered}

![ overview of the Materials Project MCP-enabled architecture
integrating the user interface, LLM provider with MCP client, the
protocol layer, MCP server, and the Materials Project
API.](mcp_architecture.png){#fig:mcp_architecture width="90%"}

This MCP Host for the Materials Project serves as the application that
establishes connections with multiple MCP servers to provide users
access to distributed computational resources as in Figure
[1](#fig:mcp_architecture){reference-type="ref"
reference="fig:mcp_architecture"}. The MCP Host operates as an
orchestration layer; the MCP Host manages server lifecycle events,
handles authentication credentials, and implements protocol-compliant
communication channels between the host application and connected
servers. This architecture enables applications to leverage specialized
services through standardized protocols while abstracting the underlying
complexity from end-users.

## Tool Definition and Schema-Driven Design

This MCP server is designed focusing on asynchronous efficiency, modular
structure, and reproducibility. It relies on *aiohttp*, a non-blocking
asynchronous web framework that allows simultaneous processing when
retrieving extensive datasets or performing resource heavy materials
science computations[@odida2025aiohttp]. This project arranges the
server as a collection of modular, separately callable tools. This
design encourages the division of responsibilities, enabling the clear
arrangement of materials operations like structure retrieval, electronic
property inquiries, or thermodynamic evaluations into distinct, clearly
defined services.

### Schema-Driven Tool Definition {#schema-driven-tool-definition .unnumbered}

Every tool in the server represents a feature provided by the *mp_api*
client. Nonetheless, This MCP enhances these endpoints by utilizing
rigorous pydantic based validation for input/output schemas. Input
schemas verify parameters like *material IDs, numerical limits, and list
formats*, stopping incorrectly structured queries from reaching the
back-end. This reduces the excessive API load and guarantees prompt
error identification. Output schemas transform unrefined responses
typically containing nested arrays or intricate objects into formats
tailored for LLM interpretation and subsequent reasoning. This consists
of summaries in Markdown format, organized JSON segments when required.
This schema-first method guarantees reliable, understandable, and
uniform tool behavior and also acts as embedded documentation for the
usage and limitations of each method.

### Tool Functionality {#tool-functionality .unnumbered}

The server provides tools to assist with different research activities
in computational materials science within the project database; some
examples include:

-   Materials Discovery: search_materials allows intricate filtering
    based on chemical elements, band gap ranges, or stability
    requirements.

-   Structural Analysis: get_structure_by_id obtains comprehensive data
    on lattice geometry, symmetry, and site details using the
    identifier.

-   Electronic Attributes: get_electronic_bandstructure and
    get_electronic_dos_by_id deliver band dispersion graphs, density of
    states, and associated metadata such as bandgap measurements.

-   Thermodynamic Stability: get_thermodynamic_stability combines phase
    information to calculate stability metrics and rival phase diagrams.

-   Advanced Characteristics: Functions such as
    get_dielectric_properties, get_elastic_tensor, and
    get_magnetic_properties provide tensorial and functional material
    descriptors, arranged for scientific analysis and visualization.

### Protocol Compliance and Quality Assurance {#protocol-compliance-and-quality-assurance .unnumbered}

The server conforms closely to MCP standards, guaranteeing strong
compatibility with client agents. Important quality aspects consist of

-   Deterministic schemas: Ensure uniform I/O performance throughout
    different environments and usage sessions.

-   Environment scoping: Protects sensitive credentials (e.g.,
    MP_API_KEY) using runtime-injected variables and separate
    containers.

-   Thorough logging: Records operational details such as input values,
    execution durations, errors, and data modifications for
    repeatability.

-   User-focused error management: Provides organized and detailed
    feedback, enhancing usability for both beginners and experienced
    users.

## Implementation and Deployment

### Development Process {#development-process .unnumbered}

The development process employed an iterative and experiment-based
approach to integrate the Materials Project API and our designed MCP
server. [[Materials Project
MCP]{style="color: blue"}](https://github.com/benedictdebrah/materials-project-mcp).
Initial prototyping and evaluations were performed through Jupyter
notebooks to explore the functionality and architecture of mp_api
endpoints. This process enabled rapid testing of query parameters,
response types, and available data structures. Subsequently, specific
functions were abstracted and encapsulated using the mcp.tool decorator
from the MCP library. Every function received rigorous input validation
through Pydantic models, and the outputs were structured to conform with
MCP's defined response schema. This ensured compatibility with
downstream agents based on LLMs and clients aware of MCP. Throughout the
development process, emphasis was placed on:

-   Confirming semantic compatibility between tool inputs and API
    requirements.

-   Confirming parameter types and value ranges using `pydantic.Field`
    specifications.

-   Maintaining the clarity and expandability of separate tool
    components.

-   Unifying error management and output formatting (Markdown, JSON).

Upon implementation, the server was verified using the MCP Inspector, a
diagnostic and troubleshooting utility provided with the MCP package.
This tool was utilized to verify the registration, invocation behavior,
and schema integrity of each exposed utility before it was accessible to
any host or client environment.

### Local Execution {#local-execution .unnumbered}

In the development phase, the Materials MCP server can be run locally in
a Python environment that meets all dependency needs. Authentication is
necessary for the Materials Project API, achieved by exporting a valid
API key as the environment variable
MP_API_KEY=your-materials-project-api-key. The mcp command-line
interface is utilized to launch the server in development mode: mcp dev
server.py

### Dockerized Deployment {#dockerized-deployment .unnumbered}

For portable or production-grade deployment, a prebuilt Docker image
(benedict2002/materials-project-mcp) is available. This allows the
server to be launched in an isolated, reproducible environment,
minimizing dependency issues across different systems. The API key must
be passed at runtime using the -e environment flag: docker run
$-e \quad MP\_API\_KEY=your-api-key\ benedict2002/materials-project-mcp$.
This deployment method supports integration with cloud services, CI/CD
pipelines, or institutional computing clusters, and offers an efficient
alternative to managing local Python environments.

### Environment and Dependency Requirements {#environment-and-dependency-requirements .unnumbered}

In both deployment models, the host needs to guarantee compatibility
with the server's runtime environment, which comprises:

-   aiohttp: for performing asynchronous operations with HTTP servers

-   pydantic: for validating input/output data and enforcing schemas

-   mp_api and pymatgen: for retrieving and analyzing materials science
    information

-   mcp: for adherence to Model Context Protocol and registration of
    tools

The API key is securely confined through the MP_API_KEY environment
variable for authenticating requests to the Materials Project backend.

### Verification using MCP Inspector {#verification-using-mcp-inspector .unnumbered}

Prior to integration with MCP client, the server's configuration and
tool exposure can be validated using the MCP Inspector, a diagnostic
utility provided by the MCP development framework. Web Interface
Configuration:

1.  Open the MCP Inspector UI.

2.  Under the configuration tab, set the following parameters:

3.  Command: mcp

4.  Arguments: \[\"dev\", \"server.py\"\]

5.  Environment Variables: MP_API_KEY=your-key

6.  Click Connect to initialize the server and inspect registered tools.

## Workflow and User Interaction

This system is built in a way that does not expose a traditional CLI,
GUI, or RESTful frontend for direct user interaction. Instead,
interaction is delegated entirely to MCP-aware LLM clients. Examples
include: Claude Desktop: An LLM chat application that supports tool
invocation through MCP-compliant backends. VS Code Copilot Chat: An
IDE-based interface that allows developers to prompt LLM agents to
access domain-specific tools.

In this framework, the user formulates a natural language query (e.g.,
Show me the band structure of LiFePO$_4$), which the client translates
into a tool that get executed alongside with the MCP specification. This
involves selecting the appropriate tool, structuring the input
parameters, and receiving and parsing the output schema. The LLM acts as
an intelligent mediator, abstracting backend complexity and providing
results that are context rich results to the user.

# Evaluation {#sec:evaluation}

This section presents an evaluation of the Materials MCP server in the
context of typical use cases. All evaluations were performed with Claude
and VS code copilot, which served as the LLM client for interacting with
the server.

## Case Studies

### Case 1

For data retrieval tasks, to evaluate the Materials MCP server, we ask a
series of questions and assess the accuracy of the returned results from
the server. First, we tried to query all Stable LiFeO Compounds. The
Materials MCP server searched the database for all thermodynamically
stable compounds containing Li, Fe and O., including simple LiFeO
binary/ternary compounds, phosphate-based materials, silicates, and
other complex oxides, mixed-metal compounds.

::: tcolorbox
**Search Query:** what are all known stable compounds containing *Li,
Fe* and *O* ?\
**Search Result:** The database shows **35 stable compounds in total**.
The most technologically important is *LiFePO$_4$*, which is widely used
in lithium-ion batteries due to its excellent electrochemical properties
and thermal stability. The simple binary oxides *Li$_2$FeO$_3$*,
*LiFeO$_2$*, and *Li$_5$FeO$_4$* represent the fundamental stable phases
in the Li-Fe-O system.\
:::

In this case, we compared the traditional Materials Project API (MP-API)
with the Model Context Protocol server (MCP) for retrieving crystal
structure data. To validate that MCP provides results consistent with
MP-API, we executed the query using the MCP and summarized the outcomes
in figure above. As shown, MCP returns identical compound IDs and band
gaps, while requiring significantly less setup effort. However, MP-API
still offers more advanced filtering capabilities, such as granular
property constraints and complex query operators. These results
demonstrate that MCP is well suited for rapid, exploratory
materials-data inquiries, whereas MP-API remains preferable for highly
customized or large-scale screening tasks. The code implementation for
retrieving the stable compounds containing *Li, Fe, and O* using python
and results can be found in the supporting documents
[\[sec:supporting\]](#sec:supporting){reference-type="ref"
reference="sec:supporting"}.

::: tcolorbox
**Search Query:** are there any known ternary compounds in the Ba-Ti-O
system ?\
**Search Result:**

-   **BaTiO$_3$** (Barium Titanate): The most prominent compound, known
    for its ferroelectric and piezoelectric properties.

    -   Multiple polymorphs were retrieved, including entries such as
        `mp-5020`, `mp-5933`, and `mp-2998`.

    -   Band gaps range from **0.000 eV to 2.509 eV**, depending on the
        phase.

    -   The thermodynamically stable form (`mp-5020`) has an energy
        above hull of **0.000 eV**.

-   **Other notable Ba--Ti--O ternary phases** include:

    -   BaTi$_2$O$_5$ (`mp-3943`, `mp-555966`)

    -   BaTi$_4$O$_9$ (`mp-3175`)

    -   BaTi$_5$O$_{11}$ (`mp-27790`)

    -   BaTi$_6$O$_{12}$, BaTi$_6$O$_{13}$ (multiple structural
        variants)

    -   BaTi$_8$O$_{16}$ (stable and metastable forms, e.g.,
        `mp-2754073`)

The BaTiO system demonstrates a high degree of **structural
complexity**, particularly in titanium rich compositions, which form
extended oxide networks. Many compounds have **energy above hull values
$<$ 0.05 eV**, indicating thermodynamic stability or near-stability.

The **electronic properties** vary significantly:

-   Several phases are **metallic** (band gap $\approx$ 0.000 eV)

-   Others are **semiconducting**, with band gaps up to **2.7 eV**

-   The band gap of **BaTiO$_3$** is highly **structure-dependent**
:::

Investigating the BaTiO system through Materials Project (MP) validates
the platform's reliability in accessing a variety of chemically
important and technologically significant materials information. The
effective recognition of numerous stable and metastable phases,
comprising different polymorphs of BaTiO$_3$, showcases the server's
capability to facilitate extensive materials exploration and property
evaluation. Moreover, the variations in structure and electronics seen
in the obtained compounds emphasize the MCP server's capability to
facilitate targeted searches for materials tailored for specific
applications, like those utilized in electronic and ferroelectric
devices. In general, the MCP server demonstrates its effectiveness as a
significant resource for speeding up data-driven investigations in
intricate oxide systems. The implementation details can be found in the
supporting materials
[\[sec:supporting\]](#sec:supporting){reference-type="ref"
reference="sec:supporting"}

## Traditional API versus proposed Materials Project (MP) platform

We designed a comparative evaluation using five representative materials
science queries that cover common research workflows from basic
structural inquiries to complex property-based searches. Table
[\[tab:api_queries\]](#tab:api_queries){reference-type="ref"
reference="tab:api_queries"} shows the traditional Materials Project API
constructs alongside our proposed MP platform. The Materials MCP server
reduces the time and effort needed to access and utilize computational
materials data on the Materials Project Platform. To obtain data from
the Materials Project, one needs to configure the programming
environment, import and install necessary libraries, setting up API
keys, and dealing with complex JSON objects that needs manual parsing.
In order to alleviate these technical challenges, the Materials MCP
server deals with much of the technical work. With the Materials MCP
server, users can articulate their needs in simple language, and the
system manages the processes through input validation to stop errors.
This eliminates the trail and error cycle experienced when interfacing
with APIs.

::: adjustbox
max width=

::: tabular
p1cm p2.5cm p7.5cm p3.5cm **ID** & **Task** & **Materials Project API**
& **Proposed MCP Natural Language Query**\
**Q1** & Retrieve crystal structure &

``` {.python language="Python" basicstyle="\\ttfamily\\scriptsize"}
with MPRester(API_KEY) as mpr:
    docs = mpr.materials.summary.search(
        material_ids=["mp-149"],
        fields=["structure", "formula_pretty", "symmetry", "density", "volume"]
    )
    
    if docs:
        doc = docs[0]
        print(f"Formula: {doc.formula_pretty}")
        print(f"Symmetry: {doc.symmetry}")
        print(f"Density: {doc.density:.4f} ")
        print(f"Volume: {doc.volume:.4f} \k{u}")
        print(f"\nStructure:\n{doc.structure}")
```

& Show me the crystal structure of silicon.\
**Q2** & Get band gap of a compound &

``` {.python language="Python" basicstyle="\\ttfamily\\scriptsize"}
with MPRester(API_KEY) as mpr:
    # Search for LiFePO4 materials
    docs = mpr.materials.summary.search(
        formula="LiFePO4",
        fields=["material_id", "formula_pretty", "band_gap",]
    )
    
    print(f"\nFound {len(docs)} LiFePO4 entries in Materials Project\n")
    
    for doc in docs:
        print(f"Material ID: {doc.material_id}")
        print(f"Formula: {doc.formula_pretty}")
        print(f"Band Gap: {doc.band_gap} eV")
       
```

& What is the band gap of LiFePO$_4$?\
:::
:::

[]{#tab:api_queries label="tab:api_queries"}

Consequently, users dedicate less time to managing programming specifics
or troubleshooting errors, and more time analyzing results and applying
them to inform research choices. This implies that both specialists and
beginners can retrieve trustworthy data more quickly, transforming this
MCP server into a tool that enhances research efficiency and minimizes
obstacles in computational materials processes.

In the first case, the user retrieves the crystal structure of silicon.
Using the traditional Materials Project API, this requires importing the
MPRester client, specifying the material ID (mp-149), selecting multiple
output fields, and manually parsing the returned data. Any error in
syntax, API configuration, or field selection can interrupt the
workflow. With the Materials MCP server, the same task is completed by
simply asking, "Show me the crystal structure of silicon." The MCP
interface internally interprets the request, calls the appropriate API
endpoint, and returns the structure in a readable format. This
demonstrates how Materials Project MCP eliminates unnecessary code,
simplifies the data-retrieval process, and reduces the potential for
programming errors, making high-level materials information accessible
through natural language. The figure below shows the result of the
prompt.

<figure id="fig:silicon">
<div class="tcolorbox">
<p><strong>Structure Type:</strong> Face-centered cubic (FCC) diamond
structure</p>
<div class="tcolorbox">
<p><strong>Lattice Parameters:</strong></p>
<ul>
<li><p>All sides equal: <span
class="math inline"><em>a</em> = <em>b</em> = <em>c</em> = 3.849</span>
Å</p></li>
<li><p>All angles equal: <span
class="math inline"><em>α</em> = <em>β</em> = <em>γ</em> = 60.00<sup>∘</sup></span></p></li>
</ul>
</div>
<p><strong>Unit Cell:</strong> Contains 2 silicon atoms with reduced
formula Si</p>
<p><strong>Notes:</strong> This is the standard diamond cubic structure
characteristic of group IV semiconductors. The 60<span
class="math inline"><sup>∘</sup></span> angles indicate the primitive
rhombohedral setting of the conventional cubic cell.</p>
</div>
</figure>

The second case focuses on retrieving the band gap of $LiFePO_4$, a
widely studied cathode material for lithium-ion batteries.
Traditionally, the researcher must query the database using the compound
name, extract the relevant entry, and print the band_gap field. This not
only requires familiarity with Python and the materials project API but
also introduces room for mistakes in filtering, looping, or field
extraction. In contrast, The Materials Project MCP approach completes
the task with a single natural-language prompt: "What is the band gap of
$LiFePO_4$?" The server handles the query logic, ensures the input is
correctly structured, and returns the value directly. This highlights
the MCP server's efficiency in eliminating repetitive programming
overhead and allowing researchers to focus on interpreting material
properties rather than writing and debugging code. The figure below
shows the output of the band gap of $LiFePO_4$.

::: tcolorbox
**Band Gap:** 3.763 eV (insulator)

**Calculation Method:** Density functional theory (DFT)

::: tcolorbox
**Important Notes:**

-   DFT calculations often underestimate band gaps; experimental value
    may be larger

-   LiFePO$_4$ is an important cathode material in lithium-ion batteries

-   Wide band gap contributes to stability and safety characteristics
:::
:::

Lastly, we tested and compared the Materials MCP server to handle
compound queries, we designed an experiment comparing it directly with
the traditional Materials Project API. The full prompt used in this
evaluation is shown below.

::: tcolorbox
**Search Query:**\
Using the Materials Project data for material mp-149 (Silicon), perform
the following analyses: Retrieve and examine its electronic band
structure to confirm whether Silicon behaves as a semiconductor. Obtain
its dielectric properties and charge density information, and summarize
the key features. Compare its thermodynamic stability (e.g., energy
above hull, formation energy) with its surface properties across
relevant orientations. Assess whether Silicon is a suitable substrate
material for integrating or growing other materials, based on lattice
matching and surface stability. Provide any published references (DOIs)
from the Materials Project or related databases that discuss these
specific properties of Silicon. Return all findings in a structured
format with clear sections.
:::

This evaluation highlights how compound prompting can be handled by the
Materials MCP Server. In comparison with the traditional API calling,
which user would have to manually and independently write each query.
The Materials MCP demonstrated the ability to get all the required data
and present it to the user, and help with complex materials data
workflows. The results can be found in the Supporting Materials section.

# Summary and Conclusions {#sec:summary}

This study presents the Materials Project Model Context Protocol (MCP)
server, an innovative system that improves access to computational
materials information via natural language questions. By enclosing the
Materials Project API in an MCP-compliant framework, the server allows
large language models (LLMs) and agentic systems to engage with
materials data in a structured, user-friendly, and contextually rich
way. Our assessment shows that the MCP server effectively obtains
materials science data with accuracy similar to direct API queries,
while offering user-friendly interactions for those without programming
skills.

Through multiple case studies, the MCP server successfully addressed
intricate materials inquiries, including the identification of all
stable compounds within the Li--Fe--O and BaTiO systems. The server
exhibited a strong degree of accuracy, with all relevant descriptors
(e.g., band gap, energy above hull) aligning with direct API responses
within a slight margin. While the MCP results sometimes yielded fewer
entries because of a wider interpretation (e.g., embracing multi-element
compounds beyond ternaries), this tendency aligns with human intent and
facilitates exploratory scientific inquiry. In summary, the Materials
MCP server reduces the technical burden associated with querying
materials databases and allows researchers to focus on interpretation
rather than implementation. By enabling richer, multi-parameter
reasoning through natural-language prompts, the system contributes to
the broader accessibility of materials data, accelerates hypothesis
generation, and provides a scalable pathway for integrating AI-mediated
workflows into computational materials research.

## Challenges and Possible Solutions

Despite its promise, the Materials MCP server faces several limitations.

-   Restricted Filtering Options: In contrast to direct API requests
    that enable precise constraints (e.g., specific element counts,
    numerical ranges), the existing MCP interface does not offer
    user-adjustable filters unless they are clearly mentioned in the
    prompt. Possible solutions is to implement a schema-driven filtering
    layer that allows users to explicitly specify constraints (such as
    ranges, categories, or thresholds) within the MCP request.

-   Reliance on LLM Performance: The quality of responses is influenced
    by the LLM utilized and its ability to comprehend the server's
    schema. Inaccurate tool usage can result in incomplete or unrelated
    outcomes.

-   Latency and Throughput: While response times are typically
    satisfactory for individual queries, simultaneous access by numerous
    users can diminish performance, particularly for extensive or
    resource-heavy queries.

::: acknowledgement
We would first like to pay tribute to our dear departed colleague
Owusuaa Fobi, with whom we began this project. May her soul rest in
perfect peace. We also sincerely acknowledge the invaluable support of
Kwabena Koranteng Asiedu, and the implicit support of KNUST, all of
which contributed greatly to the realization of this work.
:::

# Funding Statement {#funding-statement .unnumbered}

This research received no specific grant from any funding agency in the
public, commercial, or not-for-profit sectors.

# Author Contributions {#author-contributions .unnumbered}

L.E.K.A. conceptualized the project, supervised the project, and
provided critical revisions on the manuscript. B.O.D. and E.P.F. carried
out the programming, data processing, and manuscript writing. All
authors reviewed and approved the final manuscript.

# Competing Interests {#competing-interests .unnumbered}

The authors declare no competing interests.

# Supporting Materials: Query Implementation Details {#sec:supporting .unnumbered}

## S.1 Manual MP API Query for Ternary Li--Fe--O Compounds {#s.1-manual-mp-api-query-for-ternary-lifeo-compounds .unnumbered}

To retrieve all thermodynamically stable ternary compounds in the
Li--Fe--O chemical system, we used the following Python code with the
`mp-api` client:

    from mp_api.client import MPRester

    API_KEY = "your-materials-project-api-key"

    with MPRester(API_KEY) as mpr:
        results = mpr.materials.summary.search(
            elements=["Li", "Fe", "O"],
            energy_above_hull=(0, 0.05)  # Stable or nearly stable
        )

    # Filter for ternary compounds (exactly 3 elements)
    filtered_results = [r for r in results if len(r.elements) == 3]

    # Extract key fields
    mp_data = [{
        "material_id": r.material_id,
        "band_gap": round(r.band_gap, 3),
        "energy_above_hull": round(r.energy_above_hull, 3),
        "formula": r.formula_pretty,
        "elements": r.elements
    } for r in filtered_results]

    print(f"Found {len(mp_data)} ternary Li-Fe-O compounds")

This query returned **37** entries containing only lithium, iron, and
oxygen.

## S.2 Natural Language Query to MCP Server {#s.2-natural-language-query-to-mcp-server .unnumbered}

The following natural language prompt was issued to the Materials MCP
server via an MCP-aware LLM client:

> *"What are the stable compounds containing lithium, iron, and
> oxygen?"*

The MCP server responded with a structured list of compounds and noted
that 35 stable entries were found. However, only 18 compounds were
explicitly shown in the LLM's response output. The returned materials
included both ternary Li--Fe--O compounds and more complex compositions
containing additional elements. A representative subset is shown below:

-   **Simple Li--Fe--O Compounds:**

    -   Li$_2$FeO$_3$ (mp-774155) -- Band gap: 0.091 eV

    -   LiFeO$_2$ (mp-757614) -- Band gap: 1.712 eV

    -   Li$_5$FeO$_4$ (mp-19511) -- Band gap: 3.084 eV

-   **Phosphate-Based Compounds:**

    -   LiFePO$_4$ (mp-19017) -- Band gap: 3.922 eV

    -   LiFeP$_2$O$_7$ (mp-19294) -- Band gap: 1.744 eV

    -   LiFe(PO$_3$)$_3$ (mp-26003) -- Band gap: 3.949 eV

    -   LiFe(PO$_3$)$_4$ (mp-25977) -- Band gap: 2.938 eV

    -   LiFe$_2$P$_3$O$_{10}$ (mp-27007) -- Band gap: 3.169 eV

-   **Silicate and Related Compounds:**

    -   Li$_2$FeSiO$_4$ (mp-18968) -- Band gap: 3.124 eV

    -   Li$_2$FeGeO$_4$ (mp-1210923) -- Band gap: 2.156 eV

-   **With Additional Transition Metals:**

    -   Li$_2$CrFeO$_4$ (mp-756242) -- Band gap: 2.038 eV

    -   Li$_4$Cr$_3$FeO$_8$ (mp-755721) -- Band gap: 2.354 eV

    -   LiNbFeO$_4$ (mp-756415) -- Band gap: 2.242 eV

    -   Li$_5$Ti$_6$FeO$_{16}$ (mp-769637) -- Band gap: 2.703 eV

-   **Other Complex Oxides:**

    -   LiFeBO$_3$ (mp-777049) -- Band gap: 3.162 eV

    -   LiFeAsO$_4$ (mp-510641) -- Band gap: 2.286 eV

    -   LiFeSbO$_4$ (mp-774920) -- Band gap: 1.355 eV

    -   LiFe(CO$_3$)$_2$ (mp-753899) -- Band gap: 2.103 eV

# Compound Query

## Silicon: Electronic, Dielectric, and Thermodynamic Properties

**Electronic Structure**

-   Silicon exhibits characteristic semiconductor behavior with a
    calculated band gap.

-   The electronic density of states (DOS) confirms its semiconductor
    nature, showing a clear band gap.

-   Charge density analysis reveals polarization, indicating a
    non-uniform electron distribution.

**Dielectric Properties**

-   Total dielectric constant: $\varepsilon = 13.00$ (isotropic tensor).

-   The dielectric response is fully electronic with no ionic
    contribution.

-   Refractive index: $n = 3.61$, consistent with experimental values
    for silicon.

**Thermodynamic Stability and Surface Properties**

-   Energy above hull: $0.0~\text{eV/atom}$, indicating perfect
    thermodynamic stability.

-   Multiple stable surface terminations were identified:

    -   \(100\) surface: $1.28~\text{J/m}^2$ (reconstructed).

    -   \(111\) surface: $1.30~\text{J/m}^2$ (reconstructed).

-   The most stable surfaces exhibit reconstruction to enhance
    stability.

**Substrate Suitability** Silicon demonstrates strong potential as a
substrate material due to excellent lattice matching with several
semiconductors (e.g., GaAs, Ge). The most compatible substrates include:

-   CeO$_2$: lowest mismatch energy, $1.37\times10^{-7}$.

-   BaF$_2$: low mismatch energy, $6.73\times10^{-5}$.

-   NdGaO$_3$: mismatch energy of $1.97\times10^{-5}$.

These properties make Si particularly suitable for (100) and (111)
orientations.

**Reference Information**

-   DOI: `10.17188/1190959`

-   Published by: Kristin Persson (2014)

-   Calculations performed using density functional theory (DFT)

-   Part of the Materials Project database

**Key Applications**

-   Semiconductor devices (due to its ideal band gap)

-   Substrate for epitaxial growth (due to strong lattice matching
    capabilities)

-   Optical applications (due to favorable dielectric response)

-   Electronic components (owing to structural and electronic stability)

# Data and Code Availability {#data-and-code-availability .unnumbered}

The full implementation, including source code and Docker configuration,
is available at: `github.com/benedictdebrah/materials-project-mcp`
