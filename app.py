import streamlit as st
import pandas as pd
from PIL import Image
import subprocess
import os
import base64
import pickle

# Molecular descriptor calculator
def desc_calc():
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -standardizenitro -fingerprints -descriptortypes PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, stderr=subprocess.PIPE)
    output, error = process.communicate()
    
    if error:
        st.error("Error running PaDEL-Descriptor. Ensure Java is installed and paths are correct.")
    
    if os.path.exists("molecule.smi"):
        os.remove("molecule.smi")

# File download function
def filedownload(df):
    csv = df.to_csv(index=False)
    b64 = base64.b64encode(csv.encode()).decode()
    href = f'<a href="data:file/csv;base64,{b64}" download="prediction.csv">Download Predictions</a>'
    return href

# Model building function
def build_model(input_data):
    try:
        with open('keap1_model.pkl', 'rb') as file:
            load_model = pickle.load(file)
        
        prediction = load_model.predict(input_data)
        st.header('**Prediction output**')
        
        prediction_output = pd.Series(prediction, name='pIC50')
        molecule_name = pd.Series(load_data.iloc[:, 0], name='molecule_name')
        df = pd.concat([molecule_name, prediction_output], axis=1)
        
        st.write(df)
        st.markdown(filedownload(df), unsafe_allow_html=True)
    except Exception as e:
        st.error(f"Error in model prediction: {e}")

# Streamlit UI
st.markdown("# Bioactivity Prediction App (Keap-1 Inhibitor)")
st.sidebar.header('1. Upload your CSV data')

uploaded_file = st.sidebar.file_uploader("Upload your input file", type=['txt'])

if st.sidebar.button('Predict'):
    if uploaded_file is not None:
        load_data = pd.read_table(uploaded_file, sep=' ', header=None)
        load_data.to_csv('molecule.smi', sep='\t', header=False, index=False)
        
        st.header('**Original input data**')
        st.write(load_data)

        with st.spinner("Calculating descriptors..."):
            desc_calc()

        if os.path.exists('descriptors_output.csv'):
            st.header('**Calculated molecular descriptors**')
            desc = pd.read_csv('descriptors_output.csv')
            st.write(desc)
            st.write(desc.shape)

            st.header('**Subset of descriptors from previously built models**')
            
            if os.path.exists('descriptor_list.csv'):
                Xlist = list(pd.read_csv('descriptor_list.csv').columns)
                desc_subset = desc[Xlist]
                st.write(desc_subset)
                st.write(desc_subset.shape)
                
                build_model(desc_subset)
            else:
                st.error("Missing 'descriptor_list.csv'. Ensure it is in the correct directory.")
        else:
            st.error("PaDEL-Descriptor failed. Check if Java is installed.")
    else:
        st.error("Please upload a valid text file.")
else:
    st.info('Upload input data in the sidebar to start!')
