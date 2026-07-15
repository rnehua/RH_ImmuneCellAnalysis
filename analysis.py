## Teiko Technical (Part 2: Analysis portion)
# Rene Huerta

# Importing packages needed for analysis and visualization
from pathlib import Path
from random import sample
import sqlite3
import pandas as pd
from scipy.stats import mannwhitneyu
from statsmodels.stats.multitest import multipletests

# Identifying the folder that has this script to allow program to work regardless of where repository is cloned
rootDir = Path(__file__).parent
# Defining path to the SQLite database created by load_data.py
dbPath = rootDir / "cell-count.db"

# Function for connecting to the SQLite database
def connectDB():
    # returns active connection for use in analysis and other queries
    return sqlite3.connect(dbPath)

# Part 2:Creating summary table

# Function to calcculate relative frequency of each immune cell population by sample
def freqSummary(conn):
    # Using expressions to reshape cell count data and calculate sample's total count
    query = """
        WITH population_counts AS (
            SELECT
                sample_id AS sample,
                'b_cell' AS population,
                b_cell AS count
            FROM cellCounts

            UNION ALL

            SELECT
                sample_id AS sample,
                'cd8_t_cell' AS population,
                cd8_t_cell AS count
            FROM cellCounts

            UNION ALL

            SELECT
                sample_id AS sample,
                'cd4_t_cell' AS population,
                cd4_t_cell AS count
            FROM cellCounts

            UNION ALL

            SELECT
                sample_id AS sample,
                'nk_cell' AS population,
                nk_cell AS count
            FROM cellCounts

            UNION ALL

            SELECT
                sample_id AS sample,
                'monocyte' AS population,
                monocyte AS count
            FROM cellCounts
        ),

         -- Summing the five population counts to get total count for each sample

        sample_totals AS (
            SELECT
                sample_id AS sample,
                b_cell + cd8_t_cell + cd4_t_cell + nk_cell + monocyte AS total_count
            FROM cellCounts
        )

        SELECT
            population_counts.sample,
            sample_totals.total_count,
            population_counts.population,
            population_counts.count,
       
         -- Joining population counts to matching sample totals to 
         -- calculate relative frequency of each population in each sample

            ROUND(
                100.0 * population_counts.count
                / sample_totals.total_count, 2) AS percentage

        FROM population_counts

        INNER JOIN sample_totals
            ON population_counts.sample = sample_totals.sample

        ORDER BY
            population_counts.sample,
            population_counts.population
    """
    # Executing SQL query and returning results as a pandas DataFrame
    frequencySummary = pd.read_sql_query(query, conn)

    return frequencySummary



# Part 3: Statistical Analysis

# Function for selecting melanoma PBMC samples from patients treated
# with miraclib and comparing responders with non-responders
def treatmentAnalysis(conn):

    query = """
        WITH population_frequencies AS (

            SELECT
                Samples.sample_id AS sample,
                Samples.subject_id AS subject,
                Subjects.response,
                'b_cell' AS population,
                cellCounts.b_cell AS count,

                100.0 * cellCounts.b_cell
                / (
                    cellCounts.b_cell
                    + cellCounts.cd8_t_cell
                    + cellCounts.cd4_t_cell
                    + cellCounts.nk_cell
                    + cellCounts.monocyte
                ) AS percentage

            FROM Samples

            INNER JOIN Subjects
                ON Samples.subject_id = Subjects.subject_id

            INNER JOIN cellCounts
                ON Samples.sample_id = cellCounts.sample_id

            WHERE LOWER(Subjects.condition) = 'melanoma'
              AND LOWER(Subjects.treatment) = 'miraclib'
              AND UPPER(Samples.sample_type) = 'PBMC'
              AND LOWER(Subjects.response) IN ('yes', 'no')

            UNION ALL

            SELECT
                Samples.sample_id,
                Samples.subject_id,
                Subjects.response,
                'cd8_t_cell',
                cellCounts.cd8_t_cell,

                100.0 * cellCounts.cd8_t_cell
                / (
                    cellCounts.b_cell
                    + cellCounts.cd8_t_cell
                    + cellCounts.cd4_t_cell
                    + cellCounts.nk_cell
                    + cellCounts.monocyte
                )

            FROM Samples

            INNER JOIN Subjects
                ON Samples.subject_id = Subjects.subject_id

            INNER JOIN cellCounts
                ON Samples.sample_id = cellCounts.sample_id

            WHERE LOWER(Subjects.condition) = 'melanoma'
              AND LOWER(Subjects.treatment) = 'miraclib'
              AND UPPER(Samples.sample_type) = 'PBMC'
              AND LOWER(Subjects.response) IN ('yes', 'no')

            UNION ALL

            SELECT
                Samples.sample_id,
                Samples.subject_id,
                Subjects.response,
                'cd4_t_cell',
                cellCounts.cd4_t_cell,

                100.0 * cellCounts.cd4_t_cell
                / (
                    cellCounts.b_cell
                    + cellCounts.cd8_t_cell
                    + cellCounts.cd4_t_cell
                    + cellCounts.nk_cell
                    + cellCounts.monocyte
                )

            FROM Samples

            INNER JOIN Subjects
                ON Samples.subject_id = Subjects.subject_id

            INNER JOIN cellCounts
                ON Samples.sample_id = cellCounts.sample_id

            WHERE LOWER(Subjects.condition) = 'melanoma'
              AND LOWER(Subjects.treatment) = 'miraclib'
              AND UPPER(Samples.sample_type) = 'PBMC'
              AND LOWER(Subjects.response) IN ('yes', 'no')

            UNION ALL

            SELECT
                Samples.sample_id,
                Samples.subject_id,
                Subjects.response,
                'nk_cell',
                cellCounts.nk_cell,

                100.0 * cellCounts.nk_cell
                / (
                    cellCounts.b_cell
                    + cellCounts.cd8_t_cell
                    + cellCounts.cd4_t_cell
                    + cellCounts.nk_cell
                    + cellCounts.monocyte
                )

            FROM Samples

            INNER JOIN Subjects
                ON Samples.subject_id = Subjects.subject_id

            INNER JOIN cellCounts
                ON Samples.sample_id = cellCounts.sample_id

            WHERE LOWER(Subjects.condition) = 'melanoma'
              AND LOWER(Subjects.treatment) = 'miraclib'
              AND UPPER(Samples.sample_type) = 'PBMC'
              AND LOWER(Subjects.response) IN ('yes', 'no')

            UNION ALL

            SELECT
                Samples.sample_id,
                Samples.subject_id,
                Subjects.response,
                'monocyte',
                cellCounts.monocyte,

                100.0 * cellCounts.monocyte
                / (
                    cellCounts.b_cell
                    + cellCounts.cd8_t_cell
                    + cellCounts.cd4_t_cell
                    + cellCounts.nk_cell
                    + cellCounts.monocyte
                )

            FROM Samples

            INNER JOIN Subjects
                ON Samples.subject_id = Subjects.subject_id

            INNER JOIN cellCounts
                ON Samples.sample_id = cellCounts.sample_id

            WHERE LOWER(Subjects.condition) = 'melanoma'
              AND LOWER(Subjects.treatment) = 'miraclib'
              AND UPPER(Samples.sample_type) = 'PBMC'
              AND LOWER(Subjects.response) IN ('yes', 'no')
        )

        SELECT *
        FROM population_frequencies

        ORDER BY
            population,
            response,
            sample
    """

    responseAnalysis = pd.read_sql_query(query, conn)

    return responseAnalysis

# Function to test if cell population frequencies differ
# significantly between responders and non-responders
def statisticalTests(responseAnalysis):

    # Empty list to store results for each of the immune cell populations
    results = []

    # Running a separate Mann-Whitney U test for each of the five immune cell populations
    for population in responseAnalysis["population"].unique():

        # Selecting rows belonging to the iteration's population
        populationData = responseAnalysis[
            responseAnalysis["population"] == population
        ]

        # Separating percentage values for responders/ non-responders
        responders = populationData.loc[
            populationData["response"].str.lower() == "yes",
            "percentage"
        ]

        nonResponders = populationData.loc[
            populationData["response"].str.lower() == "no",
            "percentage"
        ]

        # Using two-sided Mann-Whitney U test because we 
        # can't assume percentages are parametric
        statistic, pValue = mannwhitneyu(
            responders,
            nonResponders,
            alternative="two-sided"
        )

        # Storing results for this population in the results list
        results.append({
            "population": population,
            "responders_n": len(responders),
            "nonresponders_n": len(nonResponders),
            "responderMedian": responders.median(),
            "nonresponderMedian": nonResponders.median(),
            "uStatistic": statistic,
            "pValue": pValue
        })

    # Convert list of results into a DataFrame
    statisticalResults = pd.DataFrame(results)

    # Adjusting p-values because five separate populations are 
    # tested for more conservative test and lower the FDR
    statisticalResults["adjusted_pValue"] = multipletests(
        statisticalResults["pValue"],
        method="fdr_bh"
    )[1]

    statisticalResults["significant"] = (statisticalResults["adjusted_pValue"] < 0.05)

    return statisticalResults.sort_values("adjusted_pValue")

# Data Subset Analysis

# Function to identify baseline melanoma PBMC samples from patients who got miraclib
def baselineSamples(conn):

    query = """
        SELECT
            Projects.project_name AS project,
            Subjects.subject_id AS subject,
            Subjects.response,
            Subjects.sex,
            Samples.sample_id AS sample,
            Samples.sample_type,
            Samples.time_from_treatment_start

        FROM Samples

        INNER JOIN Subjects
            ON Samples.subject_id = Subjects.subject_id

        INNER JOIN Projects
            ON Subjects.project_id = Projects.project_id

        WHERE LOWER(Subjects.condition) = 'melanoma'
          AND LOWER(Subjects.treatment) = 'miraclib'
          AND UPPER(Samples.sample_type) = 'PBMC'
          AND Samples.time_from_treatment_start = 0

        ORDER BY
            Projects.project_name,
            Subjects.subject_id
    """

    # Running database query and storing matching baseline sample records in DataFrame
    baselineSamples = pd.read_sql_query(query, conn)

    return baselineSamples


# Function for summarizing the baseline sample characteristics
def summarizeBaselineSamples(baselineSamples):
    
    # Counting how many baseline samples came from each project
    samplesByProject = (
        baselineSamples
        .groupby("project")
        .size()
        .reset_index(name="sample_count")
    )

    # Keeping each subject once before counting responders/ non-responders
    subjectsByResponse = (baselineSamples[["subject", "response"]]
        .drop_duplicates()
        .groupby("response")
        .size()
        .reset_index(name="subject_count")
    )

    # Keeping each subject once then counting males and females
    subjectsBySex = (baselineSamples[["subject", "sex"]]
        .drop_duplicates()
        .groupby("sex")
        .size()
        .reset_index(name="subject_count")
    )
    
    # Returns three summary tables to be printed/ displayed in dashboard
    return samplesByProject, subjectsByResponse, subjectsBySex

# Function to calculate average B-cell count for male melanoma responders at baseline (t=0) 
def averageBaselineBCells(conn):

    query = """
        SELECT
            ROUND(AVG(cellCounts.b_cell), 2) AS average_b_cells

        FROM Subjects

        INNER JOIN Samples
            ON Subjects.subject_id = Samples.subject_id

        INNER JOIN cellCounts
            ON Samples.sample_id = cellCounts.sample_id

        WHERE LOWER(Subjects.condition) = 'melanoma'
          AND UPPER(Subjects.sex) = 'M'
          AND LOWER(Subjects.response) = 'yes'
          AND Samples.time_from_treatment_start = 0
    """
    result = pd.read_sql_query(query, conn)

    return result

# Main function for running all analysis portions
def main():

    # Opening database connection and then closing once queries are complete
    with connectDB() as conn:

        # Part 2: Creating + displaying the cell frequency summary
        frequencySummary = freqSummary(conn)

        print("\nPart 2: Frequency summary")
        print(frequencySummary.head(10))
        print(f"Summary shape: {frequencySummary.shape}")

        responseAnalysis = treatmentAnalysis(conn)
        statisticalResults = statisticalTests(responseAnalysis)

        # Part 3: Running statistical tests
        print("\nPart 3: Statistical results")
        print(statisticalResults)
            # CD4 T cells had an unadjusted p-value below 0.05, 
            # this association was no longer statistically significant after 
            # controlling for multiple comparisons (adjusted p = 0.0667). 
            # Therefore, we fail to reject the null hypothesis for all 
            # five immune cell populations.

        # Part 4: Summarizing baseline melanoma PBMC samples
        baselineSamplesPBMC = baselineSamples(conn)
        
        # Assigning three DataFrames returned into their own variable for printing
        (samplesByProject,subjectsByResponse, subjectsBySex) = summarizeBaselineSamples(baselineSamplesPBMC)

        print("\nPart 4: Baseline melanoma PBMC samples")
        print(baselineSamplesPBMC.head())

        print("\nSamples by project")
        print(samplesByProject)

        print("\nSubjects by response")
        print(subjectsByResponse)

        print("\nSubjects by sex")
        print(subjectsBySex)

        averageBCells = averageBaselineBCells(conn)

        print(f"Average baseline B cells of melanoma males of all sample/treatment types at baseline:\n{averageBCells}")


# Running analysis only when file directly executed
if __name__ == "__main__":
    main()