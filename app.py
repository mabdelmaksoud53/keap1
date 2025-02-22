import streamlit as st
import pandas as pd
import subprocess
import os
import joblib
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.ML.Descriptors import MoleculeDescriptors
from rdkit.Chem import Draw
from io import BytesIO
import base64

st.set_page_config(
    page_title="Virtual Screening")


st.header("Virtual Screening")
st.warning("""
[Click here to see example for input format](https://drive.google.com/file/d/1Qa-0JKVa9ccuSMd4Km6Z5nP-RLeT1Gp0/view?usp=sharing)
""")

def PUbchemfp_desc_calc():
    # Performs the descriptor calculation
    bashCommand = "java -Xms2G -Xmx2G -Djava.awt.headless=true -jar ./PaDEL-Descriptor/PaDEL-Descriptor.jar -removesalt -detectaromaticity -fingerprints  -descriptortypes ./PaDEL-Descriptor/PubchemFingerprinter.xml -dir ./ -file descriptors_output.csv"
    process = subprocess.Popen(bashCommand.split(), stdout=subprocess.PIPE, shell=True)
    output, error = process.communicate()
    os.remove('molecule.smi')
    
def descriptors(smiles):
    mols = [Chem.MolFromSmiles(s) for s in smiles] 
    calc = MoleculeDescriptors.MolecularDescriptorCalculator([x[0] for x in Descriptors._descList])
    desc_names = calc.GetDescriptorNames()
    
    Mol_descriptors =[]
    for mol in mols:

        descriptors = calc.CalcDescriptors(mol)
        Mol_descriptors.append(descriptors)
    return Mol_descriptors,desc_names 

def smiles_to_image(smiles, mol_size=(200, 200)):
    mol = Chem.MolFromSmiles(smiles)
    if mol:
        img = Draw.MolToImage(mol, size=mol_size)
        buffer = BytesIO()
        img.save(buffer, format="PNG")
        buffer.seek(0)
        img_str = base64.b64encode(buffer.read()).decode()
        return f'<img src="data:image/png;base64,{img_str}" alt="Molecule">'
    else:
        return "Invalid SMILES"	
    

# Model
def the_model(input_data):
    load_model = joblib.load('keap1_model.pkl')
    # Make prediction
    prediction = load_model.predict(input_data)
    prediction_probability=load_model.predict_proba(input_data)
    
    x=pd.DataFrame(prediction_probability,columns=["Pi","Pa"])
    st.header('Prediction Result')
    
    prediction_output = pd.Series(prediction, name='Result')

    #proba_output=pd.Series(prediction_probability,name="prediction_proba")
    
    molecule_name = pd.Series(reading_data["Molecule Name"], name='Molecule Name')
    
    Result= pd.concat([molecule_name, x, prediction_output], axis=1)
    
    result = []
    for x in Result["Result"]:
        if x==1:
            result.append("Active")
        if x==0:
            result.append("Inactive")
    Result["Result"]=result
    st.write(Result)
    prediction_csv = Result.to_csv(index=False,sep=",")
    st.download_button(label="Download prediction results",data=prediction_csv,file_name="vs_results.csv")



uplouded_file=st.file_uploader("Please upload your input file", type=['txt'])


if st.button('Predict'):
    reading_data = pd.read_table(uplouded_file, sep=' ', names=["Smiles","Molecule Name"])
    reading_data.to_csv('molecule.smi', sep = '\t', index = False, header=None)
    reading_data['Structure'] = reading_data['Smiles'].apply(smiles_to_image)
    st.subheader('Input data')
    st.markdown(reading_data.to_html(escape=False, index=False), unsafe_allow_html=True)
   ## st.write(reading_data)



    PUbchemfp_desc_calc()

 
    st.subheader('Generated PubChem_Fingerprints')
    pubfp_calc = pd.read_csv("descriptors_output.csv")
    pubfp_calc.drop('Name', axis=1, inplace=True)
    st.write(pubfp_calc)
    st.write(pubfp_calc.shape)
   
    first_column = reading_data.iloc[:, 0] 
    Mol_descriptors,desc_names =descriptors(first_column)
    df_with_200_descriptors = pd.DataFrame(Mol_descriptors,columns=desc_names)
    df=df_with_200_descriptors[["MolWt","MolLogP","NumHAcceptors","NumHDonors"]]
    st.subheader("Lipinski Rule of 5 Descriptors")
    st.write(df)
    
    
    the_model(pubfp_calc)
else:
    st.warning('Limit 250 compounds per file')
    
    