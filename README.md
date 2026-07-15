# Loblaw Bio Immune Cell Analysis

## Overview

This project analyzes immune cell population data from a clinical trial dataset. It uses SQLite for relational data storage, pandas and SciPy for analysis, Plotly for visualization, and Streamlit for the interactive dashboard.

The analysis looks at the following questions:

1. What is the relative frequency of each immune cell population in every sample?
2. Do immune cell population frequencies differ between melanoma patients who respond and do not respond to miraclib?
3. Which melanoma PBMC samples were collected at baseline from patients treated with miraclib?
4. How are those baseline samples distributed by project, treatment response, and sex?

--- 

## Repository Structure

```text
.
├── app.py
├── analysis.py
├── cell-count.csv
├── load_data.py
├── README.md
├── requirements.txt
```

### `load_data.py`

Creates the SQLite database, initializes the relational schema, loads all rows from `cell-count.csv`, and confirms the number of records inserted into each table.

### `analysis.py`

Connects to the SQLite database and performs the analyses required in Parts 2–4, including:

* calculating immune cell population frequencies,
* comparing miraclib responders and non-responders,
* running Mann–Whitney U tests,
* adjusting p-values using the Benjamini–Hochberg false discovery rate,
* identifying baseline melanoma PBMC samples,
* summarizing samples and subjects by project, response, and sex.

### `app.py`

Creates the interactive Streamlit dashboard. The dashboard contains three tabs:

* Cell Frequencies
* Response Analysis
* Baseline Subset

---

## Requirements

This project was developed using Python 3.

Install the required packages with:

```bash
pip install -r requirements.txt
```

The `requirements.txt` file should contain:

```text
pandas
scipy
statsmodels
streamlit
plotly
```
---

## Running the Project in GitHub Codespaces

### 1. Open the repository in Codespaces

From the GitHub repository page:

1. Select **Code**.
2. Select **Codespaces**.
3. Select **Create codespace on main**.

### 2. Install the dependencies

In the Codespaces terminal, run:

```bash
pip install -r requirements.txt
```

### 3. Create and load the SQLite database

Run:

```bash
python load_data.py
```

This command:

* creates `cell-count.db` in the repository root,
* creates the relational database tables,
* loads the CSV data,
* prints the number of rows inserted into each table.

### 4. Run the analysis

Run:

```bash
python analysis.py
```

This prints:

* the cell population frequency summary,
* the responder versus non-responder statistical results,
* the baseline melanoma PBMC sample subset,
* samples by project,
* subjects by response,
* subjects by sex.

### 5. Start the Streamlit dashboard

Run:

```bash
streamlit run app.py
```

In GitHub Codespaces, Streamlit will create a forwarded port. Open the forwarded URL when prompted.

If the browser does not open automatically, select the **Ports** tab in Codespaces and open the URL associated with port `8501`.

---

## Relational Database Schema

The database contains four related tables:

### Projects

| Column         | Description                       |
| -------------- | --------------------------------- |
| `project_id`   | Unique project identifier         |
| `project_name` | Unique project name               |

### Subjects

| Column       | Description                                     |
| ------------ | ----------------------------------------------- |
| `subject_id` | Unique subject identifier                       |
| `project_id` | Foreign key connecting the subject to a project |
| `condition`  | Disease or clinical indication                  |
| `age`        | Subject age                                     |
| `sex`        | Subject sex                                     |
| `treatment`  | Treatment received                              |
| `response`   | Treatment response status                       |

### Samples

| Column                      | Description                                    |
| --------------------------- | ---------------------------------------------- |
| `sample_id`                 | Unique sample identifier                       |
| `subject_id`                | Foreign key connecting the sample to a subject |
| `sample_type`               | Type of biological sample                      |
| `time_from_treatment_start` | Collection time relative to treatment start    |

### cellCounts

| Column       | Description                               |
| ------------ | ----------------------------------------- |
| `sample_id`  | Foreign key connecting counts to a sample |
| `b_cell`     | B-cell count                              |
| `cd8_t_cell` | CD8 T-cell count                          |
| `cd4_t_cell` | CD4 T-cell count                          |
| `nk_cell`    | Natural killer cell count                 |
| `monocyte`   | Monocyte count                            |

### Schema Relationships

```text
Projects
   |
   | project_id
   v
Subjects
   |
   | subject_id
   v
Samples
   |
   | sample_id
   v
cellCounts
```

---

## Schema Design Rationale

The raw CSV has project, subject, sample, and immune cell count level data in one table. The relational schema was used to separate these levels into related tables to reduce repeated data and improve effcieyncy through query organization

This design improves:

* data consistency,
* query readability,
* organization of subject and sample metadata,
* reuse of data across different analyses,
* scalability as additional projects and samples are added.

For hundreds of projects and thousands of samples, the schema can scale by adding indexes to frequently queried columns such as:

* `project_id`,
* `subject_id`,
* `sample_id`,
* `condition`,
* `treatment`,
* `response`,
* `sample_type`,
* `time_from_treatment_start`.

The database could also be extended with additional tables for:

* treatment information,
* clinical outcomes,
* longitudinal visits,
* laboratory measurements,
* genomic or molecular features,
* additional immune cell populations.

For a larger number of immune populations, the `cellCounts` table could be converted to a long-format table with the columns:

```text
sample_id
population
count
```

This would make it easier to add new population types without altering the database schema.

---

## Analysis Design

### Cell Population Frequencies

For every sample, the five cell population counts are summed to calculate the total cell count.

The relative frequency is calculated as:

```text
population count / total sample count × 100
```

The output contains:

* sample,
* total count,
* population,
* count,
* percentage.

### Responder Versus Non-Responder Analysis

The treatment-response analysis includes:

* melanoma subjects,
* patients treated with miraclib,
* PBMC samples,
* subjects with response values of `yes` or `no`.

A separate two-sided Mann–Whitney U test is run for each immune cell population.

The Mann–Whitney U test was selected because the analysis compares percentage distributions between two independent groups without requiring the relative frequencies to meet the normality assumption of a parametric two-sample t-test.

Because five separate cell populations are tested from the same dataset, the p-values are adjusted using the Benjamini–Hochberg false discovery rate procedure.

After adjustment, no immune cell population had an adjusted p-value below 0.05. Therefore, the null hypothesis was not rejected for any of the five populations.

### Baseline Subset Analysis

The baseline subset includes samples that meet all of the following conditions:

* melanoma,
* PBMC sample type,
* miraclib treatment,
* `time_from_treatment_start = 0`.

The program reports:

* the matching baseline samples,
* the number of samples from each project,
* the number of unique responder and non-responder subjects,
* the number of unique male and female subjects.

---

## Dashboard

The interactive dashboard can be launched locally with:

```bash
streamlit run app.py
```

Dashboard URL:

```text
ADD_DEPLOYED_DASHBOARD_LINK_HERE
```

If the dashboard has not yet been deployed, this placeholder should be replaced after publishing the Streamlit application.

---

## Dashboard Deployment

The dashboard can be deployed using Streamlit Community Cloud.

General steps:

1. Push the repository to GitHub.
2. Sign in to Streamlit Community Cloud.
3. Select **Create app**.
4. Choose the GitHub repository.
5. Set the main file path to:

```text
app.py
```

6. Deploy the application.
7. Copy the public URL into the Dashboard section of this README.

The deployed environment must have access to `cell-count.db`. Either commit the database file to the repository or ensure the database is created before the dashboard starts.

---

## Reproducibility Notes

* `cell-count.csv` must be located in the repository root.
* Run `load_data.py` before running `analysis.py` or `app.py`.
* The database is rebuilt each time `load_data.py` is run.
* All file paths are based on the location of the Python scripts rather than user-specific absolute paths.
* The dashboard depends on the SQLite database created by `load_data.py`.

---

## Author

Rene Huerta
