import streamlit as st
import pandas as pd
import os
import subprocess
import base64
import pickle

# Install Java (for Streamlit Cloud)
os.system("apt-get install -y default-jre")
os.system("java -version")

# Molecular descriptor calculator
def desc_calc():
    padel_command = (
        "java -Xms2G -Xmx2G -Djava.awt.headless=true "
        "-jar PaDEL-Descriptor/PaDEL-Descriptor.jar "
        "-removesalt -standardizenitro -fingerprints "
        "-descriptortypes PaDEL-Descriptor/PubchemFingerprinter.xml "
        "-dir ./ -file descriptors_output.csv"
    )
    process = subprocess.Popen(padel_command.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()

    if os.path.exists("molecule.smi"):
        os.remove("molecule.smi")

    if error:
        st.error("Error running PaDEL-Descriptor: " + error.decode())

# File download function
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building function
def build_model(input_data):
    try:
        with open("keap1_model.pkl", "rb") as model_file:
            load_model = pickle.load(model_file)
        
        prediction = load_model.predict(input_data)
        st.header("**Prediction Output**")
        
        prediction_output = pd.Series(prediction, name="pIC50")
        df = pd.concat([pd.Series(input_data.index, name="Molecule"), prediction_output], axis=1)
        
        st.write(df)
        st.markdown(filedownload(df), unsafe_allow_html=True)

    except Exception as e:
        st.error(f"Model prediction error: {str(e)}")

# Streamlit UI
st.markdown("# 🧪 Bioactivity Prediction App (Keap-1 Inhibitor)")
st.sidebar.header("1️⃣ Upload Your Molecule Data")
uploaded_file = st.sidebar.file_uploader("Upload your file (TXT format)", type=["txt"])

if uploaded_file is not None:
    try:
        load_data = pd.read_csv(uploaded_file, sep=' ', header=None)
        load_data.to_csv("molecule.smi", sep="\t", header=False, index=False)
        
        st.header("**📄 Original Input Data**")
        st.write(load_data)
        
        with st.spinner("🔬 Calculating descriptors..."):
            desc_calc()
        
        if not os.path.exists("descriptors_output.csv"):
            st.error("Descriptor calculation failed. Please check your input file.")
        else:
            st.header("**🧬 Calculated Molecular Descriptors**")
            desc = pd.read_csv("descriptors_output.csv")
            st.write(desc)
            st.write(f"Shape: {desc.shape}")
            
            if os.path.exists("descriptor_list.csv"):
                Xlist = list(pd.read_csv("descriptor_list.csv").columns)
                desc_subset = desc[Xlist]
                
                st.header("**📊 Subset of Descriptors Used for Prediction**")
                st.write(desc_subset)
                st.write(f"Shape: {desc_subset.shape}")
                
                with st.spinner("🔄 Running Prediction..."):
                    build_model(desc_subset)
            else:
                st.error("Descriptor list file not found!")

    except Exception as e:
        st.error(f"Error processing input file: {str(e)}")
else:
    st.info("📂 Upload a TXT file in the sidebar to start!")
