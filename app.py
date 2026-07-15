## Teiko Technical (Part 4: Data Subset Analysis)
# Rene Huerta

# Importing StreamLit to build interactive web dashboard and Plotly 
# for interactive data visualizations
import streamlit as st
import plotly.express as px

# Importing database connection and analysis functions from analysis.py
from analysis import (
    connectDB,
    freqSummary,
    treatmentAnalysis,
    statisticalTests,
    baselineSamples,
    summarizeBaselineSamples
)

# Configuring browser tab title and using full screen width
st.set_page_config(
    page_title="Loblaw Bio Immune Cell Analysis",
    layout="wide"
)

# Titles and descriptions for the StreamLit web dashboard
st.title("Loblaw Bio Clinical Trial Analysis")
st.write(
    "Interactive analysis of immune cell population frequencies, "
    "miraclib treatment response, and baseline melanoma PBMC samples."
)

# Opening connection to SQLite database and running queries to retrieve data for analysis and visualization
with connectDB() as conn:
    frequencySummary = freqSummary(conn)
    responseAnalysis = treatmentAnalysis(conn)
    baselineSamples = baselineSamples(conn)

# Runnning Mann-Whitney U tests and multiple testing correction 
statisticalResults = statisticalTests(responseAnalysis)

# Assigning three DataFrames returned into their own variable
(samplesByProject, subjectsByResponse, subjectsBySex) = summarizeBaselineSamples(baselineSamples)

# Different tabs for the three parts of the analysis 
tab1, tab2, tab3 = st.tabs([
    "Cell Frequencies",
    "Response Analysis",
    "Baseline Subset"
])


with tab1:

    # Header title
    st.header("Cell Population Frequencies")

    # Dropdown list to select sample to view and filter for
    selectedSample = st.selectbox(
        "Select a sample",
        sorted(frequencySummary["sample"].unique())
    )
    selectedSampleData = frequencySummary[
        frequencySummary["sample"] == selectedSample
    ]

    # Displaying the sample's descriptive statistics
    st.dataframe(
        selectedSampleData,
        width="stretch"
    )

    # Bar chart of the selected sample's immune cell population frequencies
    frequencyChart = px.bar(
        selectedSampleData,
        x="population",
        y="percentage",
        title=f"Cell Population Frequencies for {selectedSample}",
        labels={
            "population": "Cell population",
            "percentage": "Relative frequency (%)"
        }
    )

    # Displaying the bar chart in the StreamLit dashboard
    st.plotly_chart(
        frequencyChart,
        width="stretch"
    )

# Miraclib response analysis (tab 2)

with tab2:
    # Header title
    st.header("Miraclib Response Analysis")

    # Creating boxplots to compare distributions of immune 
    # cell population frequencies between responders and non-responders
    boxplot = px.box(
        responseAnalysis,
        x="population",
        y="percentage",
        color="response",
        points="outliers",
        title=(
            "Cell Population Frequencies in Responders "
            "and Non-Responders"
        ),
        labels={
            "population": "Cell population",
            "percentage": "Relative frequency (%)",
            "response": "Response"
        }
    )

    # Displaying the boxplot in the StreamLit dashboard
    st.plotly_chart(
        boxplot,
        width="stretch"
    )

    st.subheader("Statistical Test Results")
    st.dataframe(
        statisticalResults,
        width="stretch"
    )

    # Display sample sizes, medians, Mann-Whitney U test p-values, 
    # and FDR-adjusted p-values for each immune cell population
    significantPopulations = statisticalResults[
        statisticalResults["significant"]
    ]

    # Displaying message if there's any statistically signficant cell populations 
    # with FDR-adjusted p-value below 0.05
    if significantPopulations.empty:
        st.info(
            "No cell populations had an FDR-adjusted p-value below 0.05.")
    else:
        names = ", ".join(
            significantPopulations["population"].tolist())
        st.success(
            f"Significant cell populations: {names}")

# Baseline Data Subset Analysis (tab 3)
with tab3:

    # Header title
    st.header("Baseline Melanoma PBMC Samples")

    # Displaying baseline melanoma PBMC samples from patients that got miraclib treatment
    st.dataframe(
        baselineSamples,
        width="stretch"
    )

    # Creating a bar chart showing baseline samples from each trial project
    projectChart = px.bar(
        samplesByProject,
        x="project",
        y="sample_count",
        title="Baseline Samples by Project"
    )
    st.plotly_chart(
        projectChart,
        width="stretch"
    )

    # Creating bar chart that shows the number of unique baseline subjects
    # who were responders or non-responders
    responseChart = px.bar(
        subjectsByResponse,
        x="response",
        y="subject_count",
        title="Baseline Subjects by Response"
    )
    st.plotly_chart(
        responseChart,
        width="stretch"
    )

    # Bar chart on the number of unique baseline male & female subjects in 
    # selected subset
    sexChart = px.bar(
        subjectsBySex,
        x="sex",
        y="subject_count",
        title="Baseline Subjects by Sex"
    )
    st.plotly_chart(
        sexChart,
        width="stretch"
    )
