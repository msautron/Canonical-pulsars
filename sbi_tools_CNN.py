import os
import torch
from torch.distributions import MultivariateNormal
import numpy as np
import platform
from subprocess import call,Popen,PIPE
from glob import glob
from astropy.table import Table
from astropy import units as u, constants as c
import matplotlib.pyplot as plt
from tensorflow.keras import layers, models
import tensorflow as tf
from scipy.stats import qmc,gaussian_kde
import sbi.utils as utils
from sbi.utils.user_input_checks import process_prior
from sbi.inference.base import infer
from sbi.analysis import pairplot
from sbi.utils.user_input_checks import prepare_for_sbi
from sbi.inference import simulate_for_sbi
from sbi.inference import SNPE
import re
import subprocess
from scipy import stats
from scipy.interpolate import interp1d
from astropy.coordinates import SkyCoord
import pandas as pd
from tensorflow.keras.callbacks import EarlyStopping
from tensorflow.keras.optimizers import Adam

def simulator(params):
    """
    Full pipeline code to be used by the SBI infer() function. Reads in the
    parameters for the population generation and then runs the full search and
    detect pipeline. Returns several population statistics (listed below) that 
    are used to construct the likelihood function that will be used to produce 
    posteriors.
    """
    var=''
    reg_1=re.compile("-*.{12}[|]{1}")
    X_data_PPdot_r=np.zeros((32,32,1),dtype='float32')
    X_data_PPdot_x=np.zeros((32,32,1),dtype='float32')
    #try:
    if type(params) is not dict:
        params = {
            'sigma_b'    : params[0],
            'b_mean'     : 10**(params[1]),
            'p_mean'     : params[2]*1e-3,
            'sigma_p'    : params[3],
            'BR1'    : np.round(params[4]),
            'BR2'      : np.round(params[5]),
            'BR3' : np.round(params[6]),
            'thres1' : np.round(params[7]),
            'thres2' : np.round(params[8]),
            'pdecay1': params[9],
            'pdecay2': params[10],
            't_bevol1' : 10**(params[11]),
            't_bevol2' : 10**(params[12]),
            't_bevol3' : 10**(params[13]),
            'b0_evol1' : 10**(params[14]),
            'b0_evol2' : 10**(params[15]),
            'b0_evol3' : 10**(params[16]),
            'pcst' : params[17],
            'pb' : params[18],
            'pe' : params[19],
            'A_propto' : 10**(params[20]),
            'D_propto' : 10**(params[21]),
            'M_for_K' : params[22],
            'R_for_K' : params[23]
            }

    sim_run=subprocess.run([f'bash','run_pop_for_sbi.sh', str(params['sigma_b']), str(params['b_mean']), str(params['p_mean']), str(params['sigma_p']), str(params['BR1']), str(params['BR2']), str(params['BR3']), str(params['thres1']), str(params['thres2']), str(params['pdecay1']), str(params['pdecay2']), str(params['t_bevol1']), str(params['t_bevol2']), str(params['t_bevol3']), str(params['b0_evol1']), str(params['b0_evol2']), str(params['b0_evol3']), str(params['pcst']), str(params['pb']), str(params['pe']), str(params['A_propto']), str(params['D_propto']), str(params['M_for_K']), str(params['R_for_K'])], capture_output=True,text=True)
    print("STDOUT :", sim_run.stdout)
    print("STDERR :", sim_run.stderr)
    print("Return code :", sim_run.returncode)

    with open("P_Pdot_positions.txt","r") as f:
        data=re.findall(reg_1,f.read())

    for i in range(len(data)):
        for j in range(len(data[i])-1):
            var+=data[i][j]
        data[i]=float(var)
        var=''

    P,P_dot=[],[]
    for i in range(int(len(data)/20)):
        P+=[data[20*i]]    #Period of rotation of the pulsar in seconds
        P_dot+=[data[20*i+1]]  #Period derivative of rotation of the pulsar, no units

    with open("x_file.txt","r") as f:
        data_x=re.findall(reg_1,f.read())

    for i in range(len(data_x)):
        for j in range(len(data_x[i])-1):
            var+=data_x[i][j]
        data_x[i]=float(var)
        var=''

    P_x,P_dot_x=[],[]
    for i in range(int(len(data_x)/20)):
        P_x+=[data_x[20*i]]    #Period of rotation of the pulsar in seconds
        P_dot_x+=[data_x[20*i+1]]  #Period derivative of rotation of the pulsar, no units

    #Get the density maps
    H,X,Y,_=plt.hist2d(np.log10(P),np.log10(P_dot),bins=(32,32),
                        range=((-2,np.log10(3)),(-19,-11)),cmap='RdBu')
    H=np.expand_dims(H,axis=-1)
    X_data_PPdot_r=H

    H,X,Y,_=plt.hist2d(np.log10(P_x),np.log10(P_dot_x),bins=(32,32),
                        range=((-2,np.log10(3)),(-19,-11)),cmap='RdBu')
    H=np.expand_dims(H,axis=-1)
    X_data_PPdot_x=H

    #Get the rest of the data: the scalars
    with open("info_supp.txt","r") as f:
        result1=[float(ligne.strip()) for ligne in f]

    result=np.array(result1)
    
    return result,X_data_PPdot_r,X_data_PPdot_x
    #except Exception as e:
        #print(f"Error in the simulation with parameters {params}: {e}")
        #return None,None,None
        
#Function which creates a CNN
def create_CNN():
    #Model for the NN which will find features between P-\dot{P} 
    #Density map dimension
    input_shape = (32, 32, 1)

    #Encoder -> Flatten the data and try to find feature
    encoder_input = layers.Input(shape=input_shape)
    x = layers.Conv2D(32, (3, 3), strides=1, activation='relu', padding='same')(encoder_input)
    x = layers.MaxPooling2D((2, 2), strides=2, padding='valid')(x)
    x = layers.Conv2D(64, (3, 3), strides=1, activation='relu', padding='same')(x)
    x = layers.MaxPooling2D((2, 2), strides=2, padding='valid')(x)
    x = layers.Flatten()(x)
    encoder_output = layers.Dense(32, activation='relu')(x)

    #Encoder model
    encoder = models.Model(encoder_input, encoder_output)
    #Decoder
    decoder_input = layers.Input(shape=(32,))
    x = layers.Dense(4*4*2, activation='relu')(decoder_input)
    x = layers.Reshape((4, 4, 2))(x)
    x = layers.Conv2DTranspose(32, (3, 3), strides=2, activation='relu', padding='same')(x)
    x = layers.Conv2DTranspose(16, (3, 3), strides=2, activation='relu', padding='same')(x)
    x = layers.Conv2DTranspose(8, (3, 3), strides=2, activation='relu', padding='same')(x)
    decoder_output = layers.Conv2D(1, (3, 3), activation='sigmoid', padding='same')(x)

    #Decoder model
    decoder = models.Model(decoder_input, decoder_output)

    # Autoencoder (encoder + decoder)
    autoencoder_input = encoder_input
    autoencoder_output = decoder(encoder_output)
    autoencoder = models.Model(autoencoder_input, autoencoder_output)

    #Compile the autoencoder model
    #autoencoder.compile(optimizer='adam', loss='mse')
    autoencoder.compile(Adam(learning_rate=1e-4), loss='mse')
    return encoder,decoder,autoencoder
    
#Run the simulation for some time to obtain data for SBI
def Repeat_sim_and_save(num_sim,num_training,prior) :
    """
    num_sim : Choose the number of simulation for the density estimator
    num_training : Choose the number of simulation for the training of the CNN
    prior : provide the prior 
    This function allows to save the data for a chosen number of simulations in order
    to use SBI in the future.
    """
    #Deactivate GPU utilization for the python code
    tf.config.set_visible_devices([], 'GPU')

    #LHS for inference
    sampler = qmc.LatinHypercube(d=24)
    sample_params = sampler.random(n=num_sim)
    l_bounds=np.array(prior.support.base_constraint.lower_bound)
    u_bounds=np.array(prior.support.base_constraint.upper_bound)
    sample_scaled = qmc.scale(sample_params,l_bounds,u_bounds)
    sample_scaled = sample_scaled.astype(np.float32)
    	
    #LHS for training
    sample_params = sampler.random(n=num_training)
    sample_scaled_training = qmc.scale(sample_params,l_bounds,u_bounds)
    sample_scaled_training = sample_scaled_training.astype(np.float32)
    	
    #Training of the CNN
    nb_removed=0
    for k in range(num_training):
        params=sample_scaled_training[k-nb_removed,:]
        tab1,tab2,tab3=np.zeros((1,5),dtype='float32'), \
                                        np.zeros((32,32,1),dtype='float32'), \
                                        np.zeros((32,32,1),dtype='float32')
        tab1,tab2,tab3=simulator(params)
        if tab1 is not None:
                with open('result_training.txt', 'a+') as f:
                    for element in tab2:
                        for element2 in element:
                            for element3 in element2:
                                f.write(f'{element3} ')
                    for element in tab3:
                        for element2 in element:
                            for element3 in element2:
                                f.write(f'{element3} ')
                    f.write('\n')
                with open('params_training.txt','a+') as f:
                    for element in params:
                        f.write(f'{element} ')
                    f.write('\n')
        else : 
            nb_removed+=1
            #Remove the parameters that produced no results
            sample_scaled_training = np.delete(sample_scaled_training, 
                                                k-nb_removed, axis=0)

    #Simulate the normal pulsar population with the input parameters sampled with the LHS method
    nb_removed=0
    for k in range(num_sim):
        params=sample_scaled[k-nb_removed,:]
        tab1,tab2,tab3=np.zeros((1,5),dtype='float32'), \
                                      np.zeros((32,32,1),dtype='float32'),\
                                      np.zeros((32,32,1),dtype='float32')
        tab1,tab2,tab3=simulator(params)
        if tab1 is not None:
            with open('result_inference.txt','a+') as f:
                for element in tab1:
                    f.write(f'{element} ')
                for element in tab2:
                    for element2 in element:
                        for element3 in element2:
                            f.write(f'{element3} ')
                for element in tab3:
                    for element2 in element:
                        for element3 in element2:
                            f.write(f'{element3} ')
                f.write('\n')
            with open('params_inference.txt','a+') as f:
                for element in params:
                    f.write(f'{element} ')
                f.write('\n')
        else : 
            nb_removed+=1
            #Remove the parameters that produced no results
            sample_scaled = np.delete(sample_scaled, k-nb_removed, axis=0) 
    
    print(f'Number of training simulations: {num_training}\n Number of simulations for the inference: {num_sim}')
    print('----------DONE----------')

#Run the SBI pipeline thanks to the data stored in result_inference and result_training
def SBI_from_datafiles(prior,validation,show):
    '''
    This function allows to do the inference with data stored in txt files for the training of the CNN and the inference part
    It returns the posterior obtained on the input parameters
    prior: Prior used for the data
    validation: boolean indicating if we are validating the pipeline (choose True) or infering parameters by comparing it to the observations (choose False)
    show: boolean indicating if you want to have the corner plot shown in the end (choose True)
    '''
    #Deactivate parallelization on python
    tf.config.set_visible_devices([], 'GPU')

    #Get the results from the data files
    num_train=0
    with open("result_training.txt","r") as f:
        lines = f.readlines()
        f.seek(0)

    num_train=len(lines)

    X_train_PPdot_r=np.zeros((num_train,32,32,1))
    X_train_PPdot_x=np.zeros((num_train,32,32,1))

    all_matrix=[]
    with open('result_training.txt', 'r') as f:
        for line in f:
            data = np.array(list(map(float, line.split())))
            matrix_line = []
            for i in range(2):
                start = i * 32 * 32
                end = (i + 1) * 32 * 32
                matrix = data[start:end].reshape((32, 32, 1))
                matrix_line.append(matrix)
            all_matrix.append(matrix_line)

    for i in range(num_train):
        X_train_PPdot_r[i]=all_matrix[i][0]
        X_train_PPdot_x[i]=all_matrix[i][1]

    num_infer=0
    with open("result_inference.txt","r") as f:
        lines = f.readlines()
        f.seek(0)

    num_infer=len(lines)

    X_data_PPdot_r=np.zeros((num_infer,32,32,1))
    X_data_PPdot_x=np.zeros((num_infer,32,32,1))
    result=np.zeros((num_infer,69),dtype='float32')

    all_matrix,all_matrix_float=[],[]
    with open('result_inference.txt', 'r') as f:
        for line in f:
            data = np.array(list(map(float, line.split())))
            matrix_1_5 = data[:5].reshape((1, 5))
            remaining_data = data[5:]
            all_matrix_float.append(matrix_1_5)
            matrix_line=[]
            for i in range(2):
                start = i * 32 * 32
                end = (i + 1) * 32 * 32
                matrix = remaining_data[start:end].reshape((32, 32, 1))
                matrix_line.append(matrix)
            all_matrix.append(matrix_line)

    for i in range(num_infer):
        result[i,0:5]=all_matrix_float[i]
        X_data_PPdot_r[i]=all_matrix[i][0]
        X_data_PPdot_x[i]=all_matrix[i][1]

    with open("params_training.txt","r") as f:
        lines=f.readlines()

    params_training=np.zeros((num_train,24))
    for i,line in enumerate(lines):
        values=line.strip().split()
        params_training[i,:]=np.array(values,dtype=float)

    with open("params_inference.txt","r") as f:
        lines=f.readlines()

    params_inference=np.zeros((num_infer,24))
    for i,line in enumerate(lines):
        values=line.strip().split()
        params_inference[i,:]=np.array(values,dtype=float)

    sample_train_torch=torch.tensor(params_training,dtype=torch.float32)
    sample_infer_torch=torch.tensor(params_inference,dtype=torch.float32)

    #Train the CNNs
    encoder_P_Pdot_r,decoder_P_Pdot_r,autoencoder_P_Pdot_r=create_CNN()
    encoder_P_Pdot_x,decoder_P_Pdot_x,autoencoder_P_Pdot_x=create_CNN()

    early_stop = EarlyStopping(
    monitor='val_loss',
    patience=20,
    min_delta=1e-4,
    restore_best_weights=True
    )

    autoencoder_P_Pdot_r.fit(X_train_PPdot_r, X_train_PPdot_r, epochs=50, batch_size=25, validation_split=0.2)
    autoencoder_P_Pdot_x.fit(X_train_PPdot_x, X_train_PPdot_x, epochs=50, batch_size=25, validation_split=0.2)

    features_result_PPdot_r=encoder_P_Pdot_r.predict(X_data_PPdot_r)
    features_result_PPdot_x=encoder_P_Pdot_x.predict(X_data_PPdot_x)
    for i in range(len(features_result_PPdot_r)):
        result[i,5:37] = features_result_PPdot_r[i,:]
        result[i,37:69]= features_result_PPdot_x[i,:]

    result_infer_torch=torch.tensor(result,dtype=torch.float32)

    #Train the inference density estimator
    inference = SNPE(prior=prior)
    density_estimator = inference.append_simulations(sample_infer_torch, result_infer_torch).train()

    #Get posterior
    posterior = inference.build_posterior(density_estimator)

    #Validation part
    if (validation==True):
        array_validation_5params,X_valid_PPdot_r,X_valid_PPdot_x=simulator(np.array([0.5,np.log10(2.75e8),129,0.45,25,45,80,8000,50000,0.3,0.58,np.log10(2e5),np.log10(5e4),np.log10(7e4),np.log10(2e9),np.log10(1e8),np.log10(3e8),26.15,0.06,0.6,np.log10(1.47e9),np.log10(883.1),1.4,12000]))
        #Get the density maps from the observations
        density_map_matPPdot_r=np.zeros((num_infer,32,32,1))
        for i in range(num_infer):
            density_map_matPPdot_r[i,:,:,:]=X_valid_PPdot_r
        features_PPdot_r=encoder_P_Pdot_r.predict(density_map_matPPdot_r)

        print(features_PPdot_r[0])

        reconstructed_map = decoder_P_Pdot_r.predict(features_PPdot_r)
        print(X_valid_PPdot_r.shape)
        print(reconstructed_map.shape)
        plt.figure(figsize=(10,4))
        plt.subplot(1,2,1)
        plt.title("Original")
        plt.imshow(X_valid_PPdot_r[:,:,0])
        plt.subplot(1,2,2)
        plt.title("Reconstruction")
        plt.imshow(reconstructed_map[0,:,:,0])
        plt.savefig('CNN_val_r.pdf',dpi=300)

        density_map_matPPdot_x=np.zeros((num_infer,32,32,1))
        for i in range(num_infer):
            density_map_matPPdot_x[i,:,:,:]=X_valid_PPdot_x
        features_PPdot_x=encoder_P_Pdot_x.predict(density_map_matPPdot_x)

        print(features_PPdot_x[0])

        reconstructed_map_x = decoder_P_Pdot_x.predict(features_PPdot_x)
        print(X_valid_PPdot_x.shape)
        print(reconstructed_map_x.shape)
        plt.figure(figsize=(10,4))
        plt.subplot(1,2,1)
        plt.title("Original")
        plt.imshow(X_valid_PPdot_x[:,:,0])
        plt.subplot(1,2,2)
        plt.title("Reconstruction")
        plt.imshow(reconstructed_map_x[0,:,:,0])
        plt.savefig('CNN_val_x.pdf',dpi=300)

        # "make an observation" of known population
        observation = torch.cat((torch.tensor(array_validation_5params), torch.tensor(features_PPdot_r[0, :]),torch.tensor(features_PPdot_x[0,:])), dim=0)

        # sample the posterior
        samples = posterior.sample((100000,),x=observation)
        median=np.median(samples,axis=0)
        tolerance_more,tolerance_less,best_estimate=[],[],[]
        samples_numpy=samples.numpy()
        for i in range(samples_numpy.shape[1]):
            column_data = samples_numpy[:, i]

            # Creation of an histogram to approximate the distribution
            counts, bin_edges = np.histogram(column_data, bins=100, density=True)

            # Mode computation
            max_index = np.argmax(counts)
            mode_val_abscisse = (bin_edges[max_index] + bin_edges[max_index + 1]) / 2
            mode_val_ordonnee = counts[max_index]
            best_estimate.append(mode_val_abscisse)

            #Confidence interval at 95% computation
            tolerance_less.append(np.percentile(column_data, 2.5))  # 2.5 % for - limit
            tolerance_more.append(np.percentile(column_data, 97.5))  # 97.5 % for + limit 

        with open("best_estimate.txt","w+") as f:
            f.write(" ".join(map(str, best_estimate)) + "\n")
            f.write(" ".join(map(str, tolerance_less)) + "\n")
            f.write(" ".join(map(str, tolerance_more)) + "\n")
        params_known=np.array([0.5,np.log10(2.75e8),129,0.45,25,45,80,8000,50000,0.3,0.58,np.log10(2e5),np.log10(5e4),np.log10(7e4),np.log10(2e9),np.log10(1e8),np.log10(3e8),26.15,0.06,0.6,np.log10(1.47e9),np.log10(883.1),1.4,12000])

    #Comparison with observations part
    elif (validation==False):
        #Reading in the observations of normal pulsars for comparison to simulated data
        #Getting the scalar info
        with open("info_supp_obs.txt","r") as f:
            scalars_obs=[float(ligne.strip()) for ligne in f]

        temp=scalars_obs[1]
        scalars_obs[1]=scalars_obs[3]
        scalars_obs[3]=temp

        print(len(scalars_obs))
        print(scalars_obs)

        #Getting the info about the P-Pdot for radio pulsars
        reg_3=re.compile("[-+]?\d*[.]\d*[Ee]*[-+]*\d*")
        reg_survey=re.compile(r"\b(?!NULL\b)[\w,]{3,}\b")

        with open("fast_fermi_pmps.txt","r") as f:
            data2=re.findall(reg_3,f.read())

        with open("fast_fermi_pmps.txt","r") as f:
            data_survey=re.findall(reg_survey,f.read())

        type_pulsar_obs=[]
        for i in range(len(data_survey)):
            if ("fast" in data_survey[i] or "pks" in data_survey[i]) and not "Fermi" in data_survey[i]:
                type_pulsar_obs.append(1) #Pulsar radio
            elif not ("fast" in data_survey[i] and not "pks" in data_survey[i]) and "Fermi" in data_survey[i]:
                type_pulsar_obs.append(2) #Pulsar gamma
            elif ("fast" in data_survey[i] or "pks" in data_survey[i]) and "Fermi" in data_survey[i]:
                type_pulsar_obs.append(3) #Pulsar rg

        Pa,P_dota,da,E_dota=[],[],[],[]
        for i in range(int(len(data2)/9)):
            Pa+=[float(data2[i*9])] #Period of the rotation of the pulsar in seconds
            P_dota+=[float(data2[i*9+1])] #Period derivative of the rotation of the pulsar no units
            da+=[float(data2[i*9+2])] #Distance to us in kpc
            E_dota+=[float(data2[i*9+8])] #Spin down power of the pulsar in ergs/s

        P2,P_dot2,type_pulsar_obs2=[],[],[]
        for i in range(len(P_dota)):
            if P_dota[i]!=0.0 and da[i]!=0.0 and E_dota[i]!=0 and da[i]<=25:
                P2+=[Pa[i]]
                P_dot2+=[P_dota[i]]
                type_pulsar_obs2.append(type_pulsar_obs[i])

        P3,P_dot3=[],[]
        for i in range(len(type_pulsar_obs2)):
            if type_pulsar_obs2[i]==1:
                P3+=[P2[i]]
                P_dot3+=[P_dot2[i]]

        #Getting the info about the P-pdot for X-ray pulsars
        df=pd.read_excel('X_ray_data_wu_et_al.ods')
        data_X=Table.from_pandas(df)
        for i in range(len(data_X['Sorting_distance'])):
            if str(data_X['Sorting_distance'][i])[0]!='<':
                data_X['Sorting_distance'][i]=float(data_X['Sorting_distance'][i])
            else:
                val=str(data_X['Sorting_distance'][i])[1:]
                data_X['Sorting_distance'][i]=val
                data_X['Sorting_distance'][i]=float(data_X['Sorting_distance'][i])
        mask=data_X['Sorting_distance'] < 25
        data_X_filtered=data_X[mask]
        data_X=data_X_filtered
        data_X['LX_upper_limit'] = [str(x).startswith('<') for x in data_X['LX']]
        data_X['LX'] = [float(str(x).replace('<','')) for x in data_X['LX']]
        data_X['LX']=data_X['LX']*1e-7

        #Get the density maps from the observations
        density_map_PPdot_r,xedges,yedges,_=plt.hist2d(np.log10(P3),np.log10(P_dot3),bins=(32,32),range=((-2,np.log10(3)),(-19,-11)),cmap='RdBu')
        density_map_PPdot_r=np.expand_dims(density_map_PPdot_r,axis=-1)
        density_map_matPPdot_r=np.zeros((num_infer,32,32,1))
        for i in range(num_infer):
            density_map_matPPdot_r[i,:,:,:]=density_map_PPdot_r
        features_PPdot_r=encoder_P_Pdot_r.predict(density_map_matPPdot_r)

        density_map_PPdot_x,xedges,yedges,_=plt.hist2d(np.log10(data_X['P']),np.log10(data_X['Pdot']),bins=(32,32),range=((-2,np.log10(3)),(-19,-11)),cmap='RdBu')
        density_map_PPdot_x=np.expand_dims(density_map_PPdot_x,axis=-1)
        density_map_matPPdot_x=np.zeros((num_infer,32,32,1))
        for i in range(num_infer):
            density_map_matPPdot_x[i,:,:,:]=density_map_PPdot_x
        features_PPdot_x=encoder_P_Pdot_x.predict(density_map_matPPdot_x)

        # "make an observation" of known population
        observation = torch.cat((torch.tensor(scalars_obs), torch.tensor(features_PPdot_r[0, :]),torch.tensor(features_PPdot_x[0,:])), dim=0)

        # sample the posterior
        samples = posterior.sample((100000,),x=observation)
        median=np.median(samples,axis=0)
        tolerance_more,tolerance_less,best_estimate=[],[],[]
        samples_numpy=samples.numpy()
        for i in range(samples_numpy.shape[1]):
            column_data = samples_numpy[:, i]

            # Creation of an histogram to approximate the distribution
            counts, bin_edges = np.histogram(column_data, bins=100, density=True)

            # Mode computation
            max_index = np.argmax(counts)
            mode_val_abscisse = (bin_edges[max_index] + bin_edges[max_index + 1]) / 2
            mode_val_ordonnee = counts[max_index]
            best_estimate.append(mode_val_abscisse)

            #Confidence interval at 95% computation
            tolerance_less.append(np.percentile(column_data, 2.5))  # 2.5 % pour la borne inférieure
            tolerance_more.append(np.percentile(column_data, 97.5))  # 97.5 % pour la borne supérieure

    with open("best_estimate.txt","w+") as f:
        f.write(" ".join(map(str, best_estimate)) + "\n")
        f.write(" ".join(map(str, tolerance_less)) + "\n")
        f.write(" ".join(map(str, tolerance_more)) + "\n")

    print(f"Dimensions de samples: {samples.shape}")
    print(f"Best estimates: {np.shape(best_estimate)}")

    #Plot
    labelss=[r'$\sigma_{\rm B}$', r"$\log(\mu_{\rm B})$ (B in T)", r"$\mu_{\rm P}$ (ms)", r'$\sigma_{\rm P}$', r'Birth spacing 1 (yr)', r'Birth spacing 2 (yr)', r'Birth spacing 3 (yr)',r'Threshold BS1',r'Threshold BS3',r'$p_{\rm decay}^1$',r'$p_{\rm decay}^2$',r'$t_{\rm bevol1}$ (yr)',r'$t_{\rm bevol2}$ (yr)',r'$t_{\rm bevol3}$ (yr)',r'$b_0^{\rm evol1}$ (T)',r'$b_0^{\rm evol2}$ (T)',r'$b_0^{\rm evol3}$ (T)',r'pcst',r'pb',r'pe',r'Cst1',r'Cst2',r'NS mass ($M_{\odot}$)',r'NS radius (m)']
    groups = [
    [0,1,2,3,4],
    [5,6,7,8,9],
    [10,11,12,13,14],
    [15,16,17,18,19],
    [20,21,22,23]
    ]
    best_estimate=np.array(best_estimate)
    log_probability = posterior.log_prob(samples, x=observation)
    for g in groups:
        pairplot(
            samples[:, g],
            fig_size=(10,8),
            upper='kde',
            diag='kde',
            points=[best_estimate[g],params_known[g]], #When SBI validation
            #points=[best_estimate[g]], #When SBI with obs
            fig_kwargs={"cmap": "viridis"},
            labels=[labelss[i] for i in g]
        )
        plt.tight_layout()
        plt.savefig(f'pairplot{g[0]}.pdf',dpi=300)
        if(show==True):
            plt.show()
        plt.close()
    #out = pairplot(samples,
    #               fig_size=(20,18),
    #               upper='kde',
    #               diag='kde',
                   #points=[best_estimate,params_known], #When doing validation 
    #               points=[best_estimate], #When using observation
    #               fig_kwargs={
                   #    "points_labels": ["Posterior best estimate"],#,"Ground truth"],
                   #    "legend": True,
                   #    "points_colors": ["yellow"]*len(best_estimate),#"green"],
                   #    "points_offdiag": {"marker": "+", "markersize": 20},
                   #    "despine": {"offset": 0},
    #                    "cmap":"viridis",
    #               },
    #               labels=[r'$\sigma_{\rm B}$', r"$\log(\mu_{\rm B})$ (B in T)", r"$\mu_{\rm P}$ (ms)", r'$\sigma_{\rm P}$', r'BS1 (yr)', r'BS2 (yr)', r'BS3(yr)',r'Thres1',r'Thres2',r'pdecay1',r'pdecay2',r't_bevol1 (yr)',r't_bevol2 (yr)',r't_bevol3 (yr)',r'b0_evol1 (T)',r'b0_evol2 (T)',r'b0_evol3 (T)',r'pcst',r'pb',r'pe',r'Cst1',r'Cst2',r'NS mass ($M_{\odot}$)',r'NS radius (m)'])

    #plt.tight_layout()
    #plt.savefig('pairplot.pdf',dpi=300)
    #if(show==True):
    #    plt.show()
    #plt.close()
    return posterior,observation
