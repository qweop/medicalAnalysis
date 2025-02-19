import pandas as pd
from scipy.stats import chi2_contingency
import yaml


def save_histology_results(writer, results, percentages, sheet_name):
    histology_data = [
        [key, value[0], percentages[key][0], value[1], percentages[key][1], value[2], percentages[key][2]] for
        key, value in results.items()]
    df_histology = pd.DataFrame(histology_data,
                                columns=["Befund", "Gesamt (n)", "Gesamt (%)", "Männer (n)", "Männer (%)",
                                         "Frauen (n)", "Frauen (%)"])
    df_histology.to_excel(writer, sheet_name=sheet_name, index=False)


def check_multiple_types(row):
    # List of locations and gastritis types
    locations = ['A. Typ A Gastritis', 'A. Typ B Gastritis', 'A. Typ C Gastritis',
                 'K. Typ A Gastritis', 'K. Typ B Gastritis', 'K. Typ C Gastritis',
                 'D. Typ A Gastritis', 'D. Typ B Gastritis', 'D. Typ C Gastritis', ]

    # Collect unique types of gastritis present across all three locations
    types_present = set()
    for col in locations:
        if row[col] == 1:  # If gastritis of type exists in that column
            # Add the type to the set (A, B, or C) based on column name
            if 'Typ A' in col:
                types_present.add('A')
            elif 'Typ B' in col:
                types_present.add('B')
            elif 'Typ C' in col:
                types_present.add('C')

    # If there are multiple different types in at least two locations
    return len(types_present) > 1


def find_ohne_befund(df):
    # Assuming df_all is your DataFrame
    df['insgesamt_ohne_befund'] = (
        ((df['A. Ohne Befund'] == 1) &
         (df['K. Ohne Befund'] == 1) &
         (df['D. Ohne Befund'] == 1)).astype(int)
    )


def find_multiple_gastritis(df, biopsy_sites2):
    # Create a new column to count gastritis occurrences for each person
    df['gastritis_count_all'] = df[biopsy_sites2['Duodenum'] + biopsy_sites2['Antrum'] + biopsy_sites2['Corpus']].sum(
        axis=1)
    df['gastritis_count_duodenum'] = df[biopsy_sites2['Duodenum']].sum(axis=1)
    df['gastritis_count_antrum'] = df[biopsy_sites2['Antrum']].sum(axis=1)
    df['gastritis_count_corpus'] = df[biopsy_sites2['Corpus']].sum(axis=1)

    # Filter out individuals with multiple gastritis diagnoses
    df_multiple_gastritis_duodenum = df[df['gastritis_count_duodenum'] > 1]
    df_multiple_gastritis_antrum = df[df['gastritis_count_antrum'] > 1]
    df_multiple_gastritis_corpus = df[df['gastritis_count_corpus'] > 1]
    df_multiple_gastritis = df[df['gastritis_count_all'] > 1]

    return df_multiple_gastritis, df_multiple_gastritis_duodenum, df_multiple_gastritis_antrum, df_multiple_gastritis_corpus


def histology_result_counts(df, df_male, df_female, site):
    site_results = {}
    biopsy_sites = {
        'Duodenum': ['D. Ohne Befund', 'D. Intraepiteliale Lymphozytose', 'D. Peptischer Schleimhautschaden',
                     'D. Tubuläres Adenom'],
        'Antrum': ['A. Ohne Befund', 'A. Chronische Gastritis', 'A. Typ A Gastritis', 'A. Typ B Gastritis',
                   'A. Typ C Gastritis', 'A. Intestinale Metaplasie'],
        'Corpus': ['K. Ohne Befund', 'K. Chronische Gastritis', 'K. Typ A Gastritis', 'K. Typ B Gastritis',
                   'K. Typ C Gastritis', 'K. Intestinale Metaplasie']
    }
    for condition in biopsy_sites[site]:
        site_results[condition] = [
            df[condition].sum(),
            df_male[condition].sum(),
            df_female[condition].sum()
        ]

    total_counts = [len(df), len(df_male), len(df_female)]

    site_results_percent = {key: [(val[i] / total_counts[i]) * 100 if total_counts[i] > 0 else 0 for i in range(3)]
                            for key, val in site_results.items()}

    return site_results, site_results_percent


def calculate_contingency_tables_by_gender(df_all, df_male, df_female, column1, column2):
    return {
        "Overall": pd.crosstab(df_all[column1], df_all[column2]),
        "Male": pd.crosstab(df_male[column1], df_male[column2]),
        "Female": pd.crosstab(df_female[column1], df_female[column2])
    }


def calculate_chi2_correlation(df1, df2):
    contingency_table = pd.crosstab(df1, df2)
    chi2, p_value, _, _ = chi2_contingency(contingency_table)
    return p_value


def chi2_by_gender(df_all, df_male, df_female, column1, column2):
    return {
        "Overall": calculate_chi2_correlation(df_all[column1], df_all[column2]),
        "Male": calculate_chi2_correlation(df_male[column1], df_male[column2]),
        "Female": calculate_chi2_correlation(df_female[column1], df_female[column2])
    }


"""
Load config file
"""
with open("config.yaml", "r", encoding="utf-8") as file:
    config = yaml.safe_load(file)

FILE_PATH = config["file_path"]
SHEET_NAME = config["sheet_name"]
BIOPSY_SITES = config["biopsy_sites"]
LOCATION_GASTRITIS_COMBINATIONS = config["location_gastritis_combinations"]

"""
Load and modify Dataset
"""
df_all = pd.read_excel(FILE_PATH, sheet_name=SHEET_NAME, engine="openpyxl").iloc[:-1]
df_all.columns = df_all.columns.str.strip()
find_ohne_befund(df_all)
df_male = df_all[df_all["Gender"] == 2]
df_female = df_all[df_all["Gender"] == 1]

"""
General Data
"""
mean_age = df_all['Age (y)'].mean()
mean_age_male = df_male['Age (y)'].mean()
mean_age_female = df_female['Age (y)'].mean()
total_all = df_all.shape[0]
total_male = df_male.shape[0]
total_female = df_female.shape[0]
gender_distribution = df_all['Gender'].value_counts(normalize=True) * 100
# Total counts for all patients
total_typ_a_all = df_all[["A. Typ A Gastritis", "K. Typ A Gastritis", "D. Typ A Gastritis"]].sum().sum()
total_typ_b_all = df_all[["A. Typ B Gastritis", "K. Typ B Gastritis", "D. Typ B Gastritis"]].sum().sum()
total_typ_c_all = df_all[["A. Typ C Gastritis", "K. Typ C Gastritis", "D. Typ C Gastritis"]].sum().sum()

# Total counts for males
total_typ_a_male = df_male[["A. Typ A Gastritis", "K. Typ A Gastritis", "D. Typ A Gastritis"]].sum().sum()
total_typ_b_male = df_male[["A. Typ B Gastritis", "K. Typ B Gastritis", "D. Typ B Gastritis"]].sum().sum()
total_typ_c_male = df_male[["A. Typ C Gastritis", "K. Typ C Gastritis", "D. Typ C Gastritis"]].sum().sum()

# Total counts for females
total_typ_a_female = df_female[["A. Typ A Gastritis", "K. Typ A Gastritis", "D. Typ A Gastritis"]].sum().sum()
total_typ_b_female = df_female[["A. Typ B Gastritis", "K. Typ B Gastritis", "D. Typ B Gastritis"]].sum().sum()
total_typ_c_female = df_female[["A. Typ C Gastritis", "K. Typ C Gastritis", "D. Typ C Gastritis"]].sum().sum()

total_blutung_all = df_all[['I - Blutung']].sum().sum()
total_blutung_male = df_male[['I - Blutung']].sum().sum()
total_blutung_female = df_female[['I - Blutung']].sum().sum()

total_ulcera_all = df_all[['Ulcera']].sum().sum()
total_ulcera_male = df_male[['Ulcera']].sum().sum()
total_ulcera_female = df_female[['Ulcera']].sum().sum()

"""
Indications for Gastroscopy
"""
indications = {
    "Dyspepsie": [df_all['I - Dyspepsie'].sum(), df_male['I - Dyspepsie'].sum(), df_female['I - Dyspepsie'].sum()],
    "Reflux": [df_all['I - Reflux'].sum(), df_male['I - Reflux'].sum(), df_female['I - Reflux'].sum()],
    "Allgemeinsymptome": [df_all['I - Allgemeinsymptome'].sum(), df_male['I - Allgemeinsymptome'].sum(),
                          df_female['I - Allgemeinsymptome'].sum()],
    "Laborbefunde": [df_all['I - Laborbefunde'].sum(), df_male['I - Laborbefunde'].sum(),
                     df_female['I - Laborbefunde'].sum()],
    "Blutung": [df_all['I - Blutung'].sum(), df_male['I - Blutung'].sum(), df_female['I - Blutung'].sum()],
    "Wunsch des Patienten": [df_all['I - Wunsch des Patienten'].sum(), df_male['I - Wunsch des Patienten'].sum(),
                             df_female['I - Wunsch des Patienten'].sum()],
    "Unbekannt": [df_all['I - Unbekannt'].sum(), df_male['I - Unbekannt'].sum(), df_female['I - Unbekannt'].sum()]
}

"""
Macroscopic Results
"""
macroscopic_results = {
    "Reflux": [df_all['MR - Reflux'].sum(), df_male['MR - Reflux'].sum(), df_female['MR - Reflux'].sum()],
    "Ösophagitis/Gastritis/Entzündliche Veränderungen": [
        df_all['MR - Ösophagitis/Gastritis/Entzündliche Veränderungen'].sum(),
        df_male['MR - Ösophagitis/Gastritis/Entzündliche Veränderungen'].sum(),
        df_female['MR - Ösophagitis/Gastritis/Entzündliche Veränderungen'].sum()],
    "Ulcus / Erosionen": [df_all['MR - Ulcus / Erosionen'].sum(), df_male['MR - Ulcus / Erosionen'].sum(),
                          df_female['MR - Ulcus / Erosionen'].sum()],
    "Strukturanomalien": [df_all['MR - Strukturanomalien'].sum(), df_male['MR - Strukturanomalien'].sum(),
                          df_female['MR - Strukturanomalien'].sum()],
    "Tumor": [df_all['MR - Gutartige Raumforderungen (Polypen, Adenome, Papillome)'].sum(),
              df_male['MR - Gutartige Raumforderungen (Polypen, Adenome, Papillome)'].sum(),
              df_female['MR - Gutartige Raumforderungen (Polypen, Adenome, Papillome)'].sum()],
    "Blutung": [df_all['MR - Blutung'].sum(), df_male['MR - Blutung'].sum(), df_female['MR - Blutung'].sum()],
    "Ohne Befund": [df_all['MR - Ohne Befund'].sum(), df_male['MR - Ohne Befund'].sum(),
                    df_female['MR - Ohne Befund'].sum()],
    "Klinisch irrelevant": [df_all['MR - Klinisch irrelevant'].sum(), df_male['MR - Klinisch irrelevant'].sum(),
                            df_female['MR - Klinisch irrelevant'].sum()]
}

"""
Histological Results by Biopsy Sites
"""
duodenum_results, duodenum_percent = histology_result_counts(df_all, df_male, df_female, 'Duodenum')
antrum_results, antrum_percent = histology_result_counts(df_all, df_male, df_female, 'Antrum')
corpus_results, corpus_percent = histology_result_counts(df_all, df_male, df_female, 'Corpus')

"""
Find people with multiple gastritis conditions
"""
df_multiple_gastritis, df_multiple_gastritis_duodenum, df_multiple_gastritis_antrum, df_multiple_gastritis_corpus = \
    find_multiple_gastritis(df_all, BIOPSY_SITES)
matching_participants = df_all[df_all.apply(check_multiple_types, axis=1)]

"""
Correlations for Dyspepsia with Gastritis Types
"""
p_values_gastritis_A_A = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'A. Typ A Gastritis')
p_values_gastritis_B_A = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'A. Typ B Gastritis')
p_values_gastritis_C_A = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'A. Typ C Gastritis')

p_values_gastritis_A_K = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'K. Typ A Gastritis')
p_values_gastritis_B_K = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'K. Typ B Gastritis')
p_values_gastritis_C_K = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'K. Typ C Gastritis')

p_values_gastritis_A_D = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'D. Typ A Gastritis')
p_values_gastritis_B_D = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'D. Typ B Gastritis')
p_values_gastritis_C_D = chi2_by_gender(df_all, df_male, df_female, 'I - Dyspepsie', 'D. Typ C Gastritis')

"""
Crosstables for Participants Ohne Befund
"""
crosstable_ohne_befund_A = calculate_contingency_tables_by_gender(df_all, df_male, df_female, 'I - Dyspepsie',
                                                                  'A. Ohne Befund')
crosstable_ohne_befund_K = calculate_contingency_tables_by_gender(df_all, df_male, df_female, 'I - Dyspepsie',
                                                                  'K. Ohne Befund')
crosstable_ohne_befund_D = calculate_contingency_tables_by_gender(df_all, df_male, df_female, 'I - Dyspepsie',
                                                                  'D. Ohne Befund')
crosstable_ohne_befund_all = calculate_contingency_tables_by_gender(df_all, df_male, df_female, 'I - Dyspepsie',
                                                                    'insgesamt_ohne_befund')

"""
NSAR Use and Type C Gastritis / Blood / Ulcera correlation
"""
p_values_NSAR_C = chi2_by_gender(df_all, df_male, df_female, 'NSAR', 'A. Typ C Gastritis')
p_values_blood_NSAR = chi2_by_gender(df_all, df_male, df_female, 'I - Blutung', 'NSAR')
p_values_ulcera_NSAR = chi2_by_gender(df_all, df_male, df_female, 'Ulcera', 'NSAR')

NSAR_C_all = ((df_all["NSAR"] == 1) &
              (df_all[["A. Typ C Gastritis", "K. Typ C Gastritis", "D. Typ C Gastritis"]].sum(axis=1) >= 1)).sum()
NSAR_C_male = ((df_male["NSAR"] == 1) &
               (df_male[["A. Typ C Gastritis", "K. Typ C Gastritis", "D. Typ C Gastritis"]].sum(axis=1) >= 1)).sum()
NSAR_C_female = ((df_female["NSAR"] == 1) &
                 (df_female[["A. Typ C Gastritis", "K. Typ C Gastritis", "D. Typ C Gastritis"]].sum(axis=1) >= 1)).sum()

NSAR_Blutung_all = (df_all[["NSAR", "I - Blutung"]] == 1).all(axis=1).sum()
NSAR_Blutung_male = (df_male[["NSAR", "I - Blutung"]] == 1).all(axis=1).sum()
NSAR_Blutung_female = (df_female[["NSAR", "I - Blutung"]] == 1).all(axis=1).sum()

NSAR_Ulcera_all = (df_all[["NSAR", "Ulcera"]] == 1).all(axis=1).sum()
NSAR_Ulcera_male = (df_male[["NSAR", "Ulcera"]] == 1).all(axis=1).sum()
NSAR_Ulcera_female = (df_female[["NSAR", "Ulcera"]] == 1).all(axis=1).sum()

"""
Save To Excel
"""
with pd.ExcelWriter("output_data.xlsx", engine="openpyxl") as writer:
    """
    General Data
    """
    data = [["Mittleres Alter", mean_age, mean_age_male, mean_age_female],
            ["Participants (total)", total_all, total_male, total_female],
            ["Participants (%)", 100, gender_distribution[2], gender_distribution[1]],
            ["Gastritis A (total)", total_typ_a_all, total_typ_a_male, total_typ_a_female],
            ["Gastritis A (%)", total_typ_a_all / total_all * 100, total_typ_a_male / total_male * 100,
             total_typ_a_female / total_female * 100],
            ["Gastritis B (total)", total_typ_b_all, total_typ_b_male, total_typ_b_female],
            ["Gastritis B (%)", total_typ_b_all / total_all * 100, total_typ_b_male / total_male * 100,
             total_typ_b_female / total_female * 100],
            ["Gastritis C (total)", total_typ_c_all, total_typ_c_male, total_typ_c_female],
            ["Gastritis C (%)", total_typ_c_all / total_all * 100, total_typ_c_male / total_male * 100,
             total_typ_c_female / total_female * 100],
            ["Blutung (total)", total_blutung_all, total_blutung_male, total_blutung_female],
            ["Blutung (%)", total_blutung_all / total_all * 100, total_blutung_male / total_male * 100,
             total_blutung_female / total_female * 100],
            ["Ulcera (total)", total_ulcera_all, total_ulcera_male, total_ulcera_female],
            ["Ulcera (%)", total_ulcera_all / total_all * 100, total_ulcera_male / total_male * 100,
             total_ulcera_female / total_female * 100]]

    df_patient_data = pd.DataFrame(data, columns=["Merkmal", "Gesamt", "Männer", "Frauen"])
    df_patient_data.to_excel(writer, sheet_name="Patienten Daten", index=False)

    """
    Indications for Gastroscopy
    """
    indications_data = [
        [key, value[0], (value[0] / len(df_all)) * 100, value[1], (value[1] / len(df_male)) * 100, value[2],
         (value[2] / len(df_female)) * 100] for key, value in indications.items()]
    df_indications = pd.DataFrame(indications_data,
                                  columns=["Indikation", "Gesamt (n)", "Gesamt (%)", "Männer (n)", "Männer (%)",
                                           "Frauen (n)", "Frauen (%)"])
    df_indications.to_excel(writer, sheet_name="Indikationen", index=False)

    """
    Macroscopic Results
    """
    macroscopic_data = [
        [key, value[0], (value[0] / len(df_all)) * 100, value[1], (value[1] / len(df_male)) * 100, value[2],
         (value[2] / len(df_female)) * 100] for key, value in macroscopic_results.items()]
    df_macroscopic = pd.DataFrame(macroscopic_data,
                                  columns=["Befund", "Gesamt (n)", "Gesamt (%)", "Männer (n)", "Männer (%)",
                                           "Frauen (n)", "Frauen (%)"])
    df_macroscopic.to_excel(writer, sheet_name="Makroskopische Ergebnisse", index=False)

    """
    Histological Results by Biopsy Sites
    """
    save_histology_results(writer, duodenum_results, duodenum_percent, "Histologie Duodenum")
    save_histology_results(writer, antrum_results, antrum_percent, "Histologie Antrum")
    save_histology_results(writer, corpus_results, corpus_percent, "Histologie Corpus")

    """
    Multiple gastritis conditions
    """
    df_multiple_gastritis[['Participant'] + LOCATION_GASTRITIS_COMBINATIONS].to_excel(writer,
                                                                                      sheet_name='Multiple Gastritis',
                                                                                      index=False)

    """
    Correlations for Dyspepsia with Gastritis Types
    """
    correlation_data = []
    gastritis_types = [
        ("Typ A Gastritis", p_values_gastritis_A_A, p_values_gastritis_A_K, p_values_gastritis_A_D),
        ("Typ B Gastritis", p_values_gastritis_B_A, p_values_gastritis_B_K, p_values_gastritis_B_D),
        ("Typ C Gastritis", p_values_gastritis_C_A, p_values_gastritis_C_K, p_values_gastritis_C_D)
    ]

    for gastritis_type, p_values_A, p_values_K, p_values_D in gastritis_types:
        for prefix, p_values in [("A.", p_values_A), ("K.", p_values_K), ("D.", p_values_D)]:
            column_name = f"{prefix} {gastritis_type}"
            pos_histology_all = df_all[column_name].sum()
            pos_histology_male = df_male[column_name].sum()
            pos_histology_female = df_female[column_name].sum()
            correlation_data.append(
                [f"Dyspepsie & {column_name}", pos_histology_all, (pos_histology_all / len(df_all)) * 100,
                 p_values['Overall'],
                 pos_histology_male, (pos_histology_male / len(df_male)) * 100 if len(df_male) > 0 else 0,
                 p_values['Male'],
                 pos_histology_female, (pos_histology_female / len(df_female)) * 100 if len(df_female) > 0 else 0,
                 p_values['Female']])

    df_correlation = pd.DataFrame(correlation_data,
                                  columns=["Korrelation", "Gesamt (n)", "Gesamt (%)", "p-Wert (Gesamt)", "Männer (n)",
                                           "Männer (%)", "p-Wert (Männer)", "Frauen (n)", "Frauen (%)",
                                           "p-Wert (Frauen)"])
    df_correlation.to_excel(writer, sheet_name="Korrelationen Gastritis", index=False)

    """
    Crosstables for Participants Ohne Befund
    """
    cross_tables = {
        "A. Ohne Befund": crosstable_ohne_befund_A,
        "K. Ohne Befund": crosstable_ohne_befund_K,
        "D. Ohne Befund": crosstable_ohne_befund_D,
        "Alle ohne Befund": crosstable_ohne_befund_all,
    }

    sheet_name = "Kreuztabellen Ohne Befund"
    writer.sheets[sheet_name] = writer.book.create_sheet(sheet_name)
    start_row = 0

    for gastritis_type, tables in cross_tables.items():
        # Überschrift für die Gastritisform schreiben
        writer.sheets[sheet_name].cell(row=start_row + 1, column=1, value=f"Kreuztabellen für {gastritis_type}")

        # Position für die Gruppen (Overall, Male, Female)
        col_offset = 0  # Starten in der ersten Spalte

        for group, table in tables.items():
            # Gruppe als Überschrift schreiben
            writer.sheets[sheet_name].cell(row=start_row + 2, column=col_offset + 1, value=group)

            # Tabelle speichern
            table.to_excel(writer, sheet_name=sheet_name, startrow=start_row + 3, startcol=col_offset)

            # Platz für die nächste Gruppe
            col_offset += len(table.columns) + 2  # Abstand zwischen den Tabellen

        # Platz für die nächste Gastritisform
        start_row += max(len(table) for table in tables.values()) + 5  # Abstand zwischen den Gastritisformen

    """
    NSAR Use and Type C Gastritis / Blood / Ulcera correlation
    """
    nsar_data = [["NSAR & Typ C Gastritis", NSAR_C_all, (NSAR_C_all / len(df_all)) * 100,
                  p_values_NSAR_C['Overall'],
                  NSAR_C_male, (NSAR_C_male / len(df_male)) * 100 if len(df_male) > 0 else 0,
                  p_values_NSAR_C['Male'],
                  NSAR_C_female,
                  (NSAR_C_female / len(df_female)) * 100 if len(df_female) > 0 else 0,
                  p_values_NSAR_C['Female']],
                 ["NSAR & Blutung", NSAR_Blutung_all, (NSAR_Blutung_all / len(df_all)) * 100,
                  p_values_blood_NSAR['Overall'],
                  NSAR_Blutung_male, (NSAR_Blutung_male / len(df_male)) * 100 if len(df_male) > 0 else 0,
                  p_values_blood_NSAR['Male'],
                  NSAR_Blutung_female, (NSAR_Blutung_female / len(df_female)) * 100 if len(df_female) > 0 else 0,
                  p_values_blood_NSAR['Female']],
                 ["NSAR & Ulcera", NSAR_Ulcera_all, (NSAR_Ulcera_all / len(df_all)) * 100,
                  p_values_ulcera_NSAR['Overall'],
                  NSAR_Ulcera_male, (NSAR_Ulcera_male / len(df_male)) * 100 if len(df_male) > 0 else 0,
                  p_values_ulcera_NSAR['Male'],
                  NSAR_Ulcera_female, (NSAR_Ulcera_female / len(df_female)) * 100 if len(df_female) > 0 else 0,
                  p_values_ulcera_NSAR['Female']]
                 ]

    df_nsar = pd.DataFrame(nsar_data,
                           columns=["Korrelation", "Gesamt (n)", "Gesamt (%)", "p-Wert (Gesamt)", "Männer (n)",
                                    "Männer (%)", "p-Wert (Männer)", "Frauen (n)", "Frauen (%)", "p-Wert (Frauen)"])
    df_nsar.to_excel(writer, sheet_name="Korrelationen NSAR", index=False)
