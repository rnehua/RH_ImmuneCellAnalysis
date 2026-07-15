## Teiko Technical
# Rene Huerta

# Importing packages and functions needed for file paths, SQLite database, and CSV file reading
from pathlib import Path
import sqlite3
import pandas as pd

# Identifying the folder that has this script to allow program to work regardless of where repository is cloned
rootDir = Path(__file__).parent

# Defining paths for input CSV file and SQLite database file
csvPath = rootDir / 'cell-count.csv'
dbPath = rootDir / 'cell-count.db'

# Assigning column names of the five immune cell population from .csv file to a list for later use
cellPopulations = ["b_cell", "cd8_t_cell", "cd4_t_cell",
                    "nk_cell","monocyte",]

# Function for creating the relational database schema for the cell count data/metadata
def createSchema(conn):
    #Creating four related tables
        # Projects:
            # Table of each unique clinical trial project in the dataset
        # Subjects:
            # Table of unique subject-level information in the dataset
        # Samples:
            # Table of biological sample-level information in the dataset
        # cellCounts:
            # Table of each immune cell population count measured in a sample

    #Creating cursor object to execute SQL commands on the database and retrieving them
    cursor = conn.cursor()

    # Removing existing tables to allow the script to rebuild the database smoothly each run
    cursor.execute("DROP TABLE IF EXISTS cellCounts")
    cursor.execute("DROP TABLE IF EXISTS Samples")
    cursor.execute("DROP TABLE IF EXISTS Subjects")
    cursor.execute("DROP TABLE IF EXISTS Projects")

    # Projects table creation
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Projects (
                   project_id INTEGER PRIMARY KEY AUTOINCREMENT,
                   project_name TEXT UNIQUE NOT NULL
                   )
                """)
    
    # Subjects table creation
    # linked via project_id
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Subjects (
                subject_id TEXT PRIMARY KEY,
                project_id INTEGER NOT NULL,
                condition TEXT,
                age INTEGER,
                sex TEXT,
                treatment TEXT,
                response TEXT,
                   
                FOREIGN KEY (project_id)
                   REFERENCES Projects(project_id)
        )
    """)

    # Samples table creation
    # linked via subject_id
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS Samples (
                sample_id TEXT PRIMARY KEY,
                subject_id TEXT,
                sample_type TEXT,
                time_from_treatment_start INTEGER,

                FOREIGN KEY (subject_id)
                    REFERENCES Subjects(subject_id)
        )
    """)

    # Cell Counts table creation
    # linked via sample_id
    cursor.execute("""
        CREATE TABLE IF NOT EXISTS cellCounts (
                sample_id TEXT PRIMARY KEY,
                b_cell INTEGER,
                cd8_t_cell INTEGER,
                cd4_t_cell INTEGER,
                nk_cell INTEGER,
                monocyte INTEGER,

                FOREIGN KEY (sample_id)
                    REFERENCES Samples(sample_id)
        )
    """)    

    # Saving the schema changes to the database file        
    conn.commit()

# Function to load unique project names into Project table
def loadProjects(conn, cellCount):

    cursor = conn.cursor()

    # Retrieiving each unique project name from the CSV and removing any missing values
    projects = cellCount["project"].dropna().unique()

    # Using SQL query to insert project name into Projects table
    for project in projects:
        cursor.execute("""
            INSERT OR IGNORE INTO Projects (project_name)
            VALUES (?)
        """, (project,))

    conn.commit()

# Function to load the Subjects table
def loadSubjects(conn, cellCount):

    cursor = conn.cursor()

    # Keeping one row per subject as subjects can have multiple samples
    subjects = cellCount[
        [
            "subject",
            "project",
            "condition",
            "age",
            "sex",
            "treatment",
            "response"
        ]
    ].drop_duplicates(subset=["subject"])

    for _, row in subjects.iterrows():

        # Finding the project ID that matches the subject's project name
        cursor.execute("""
            SELECT project_id
            FROM Projects
            WHERE project_name = ?
        """, (row["project"],))

        # Retrieves first matching row as a tuple and takes integer project_id to be used as foreign key
        projectId = cursor.fetchone()[0]

        # None is converted to NULL in database for downstream use
        response = (None
            if pd.isna(row["response"])
            else row["response"]
        )

        cursor.execute("""
            INSERT INTO Subjects (
                subject_id,
                project_id,
                condition,
                age,
                sex,
                treatment,
                response
            )
            VALUES (?, ?, ?, ?, ?, ?, ?)
        """, (
            row["subject"],
            projectId,
            row["condition"],
            int(row["age"]),
            row["sex"],
            row["treatment"],
            response
        ))

    conn.commit()

# Function to load the Samples table
def loadSamples(conn, cellCount):

    cursor = conn.cursor()

    # Insert one row for each sample in the CSV into the Samples table
    for _, row in cellCount.iterrows():

        cursor.execute("""
            INSERT INTO Samples (
                sample_id,
                subject_id,
                sample_type,
                time_from_treatment_start
            )
            VALUES (?, ?, ?, ?)
        """, (
            row["sample"],
            row["subject"],
            row["sample_type"],
            int(row["time_from_treatment_start"])
        ))

    conn.commit()

# Function to load immune cell counts into the CellCounts table
def loadCellCounts(conn, cellCount):

    cursor = conn.cursor()

    # Inserting five immune cell population counts for each sample in the CSV
    for _, row in cellCount.iterrows():

        cursor.execute("""
            INSERT INTO cellCounts (
                sample_id,
                b_cell,
                cd8_t_cell,
                cd4_t_cell,
                nk_cell,
                monocyte
            )
            VALUES (?, ?, ?, ?, ?, ?)
        """, (
            row["sample"],
            int(row["b_cell"]),
            int(row["cd8_t_cell"]),
            int(row["cd4_t_cell"]),
            int(row["nk_cell"]),
            int(row["monocyte"])
        ))

    conn.commit()

# Function to confirm if rows were properly loaded into each table
def checkDB(conn):
    cursor = conn.cursor()

    tableNames = ["Projects", "Subjects", "Samples", "cellCounts"]

    # Printing number of rows in each table as extra check
    print("\nDatabase table counts:")
    for tableName in tableNames:
        cursor.execute(f"SELECT COUNT(*) FROM {tableName}")
        count = cursor.fetchone()[0]
        print(f"Number of rows in {tableName}: {count}")

# Main function for running complete database creation and loading process
def main():
    # Loading CSV into pandas DataFrame
    cellCount = pd.read_csv(csvPath)

    #Inspecting dataset to confirm it loaded correctly and reviews characteristics
    print(f"Column names: {cellCount.columns.tolist()}")
    print(f"Shape: {cellCount.shape}")
    print(f"First 5 rows:\n{cellCount.head()}")

    print("\nDataFrame information:")
    cellCount.info()    

    #Opening SQLite connection to create database file
    with sqlite3.connect(dbPath) as conn:

        createSchema(conn)
        loadProjects(conn, cellCount)
        loadSubjects(conn, cellCount)
        loadSamples(conn, cellCount)
        loadCellCounts(conn, cellCount)
        checkDB(conn)

# Running analysis only when file directly executed
if __name__ == "__main__":
    main()