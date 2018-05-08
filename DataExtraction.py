#!python3

import pandas as pd
import glob
import numpy as np
import os

# todo: Have the program detect the parameters that are in the file before anything.
#       This will remove the problems associated with files not having a 
#       homogeneous quantity and order of exported parameters.
def extraction(df, nome=''):
    df = df.replace(to_replace=r'\s+', value=np.nan, regex=True)
    
    # Curva de fluxo
    df_CF = df[['GP', 'Eta','T', 'Tau']].dropna() # Extrai CF porque os Osc tem n/a nos valores de eta
    df_CF = df_CF.add_prefix(nome)

    if df_CF.size == 0:
        df_CF = None
    
    
    df_osc = df[['w','Tau', 'G1', 'G2', 'T']].dropna()
    
    if df_osc.size != 0:  # Checa se tem dado de oscilatório
        ## Oscilatório Tensão
        contagem_w_OT = df_osc['w'].value_counts() # Conta quantas vezes um valor de freq se repetiu
        if contagem_w_OT.size == 0: # Matriz vazia
            w_mais_freq_OT = pd.Series([])
            df_OT = None # Sem qualquer dado oscilatório
        elif contagem_w_OT.max() == 1:  # Não repete nenhum valor
            df_OT = None  # Sem dado de Osc Tens
            w_mais_freq_OT = 0  # Coloca um valor para comparar
        else:
            w_mais_freq_OT = contagem_w_OT.idxmax() # Valor de freq que mais repete (1 Hz): indica Osc Tens
            indice_w_mais_freq = contagem_w_OT.max()

            df_OT_m = df_osc['w'] == w_mais_freq_OT  # Criar máscara para o Osc Tens
            df_OT = df_osc[df_OT_m]  # Aplicar a máscara.  #### Problema aqui!
            df_OT = df_OT[['Tau', 'G1', 'G2', 'T']]  # Separa só os valores de interesse

            if df_OT.index[-1] - 1 != df_OT.index[-2]: # É possível que haja uma coincidencia e repita um valor de Freq.
                #df_OT = df_OT.drop(index = df_OT.index[-1]) # Remove último ponto se ele não for contínuo com o anterior
                df_OT = df_OT.drop(df_OT.index[-1])

            df_OT = df_OT.add_prefix(nome)  # Coloca o nome para exportação

        ## Frequencia
        
        
        if contagem_w_OT.size != 0: # Revisar aqui
            df_OF_m = df_osc['w'] != w_mais_freq_OT # complementar ao Osc Tens
        else:
            df_OF_m = [] # Matriz vazia

        for i, item in enumerate(df_OF_m):
            try:
                if df_OF_m[i-1] == True and df_OF_m[i+1] == True and df_OF_m[i] == False:
                    df_OF_m[i] = True
            except KeyError:
                pass
            except IndexError:
                pass
                
        df_OF = df_osc[df_OF_m]
        df_OF = df_OF[['w', 'G1', 'G2', 'T']]
        df_OF = df_OF.add_prefix(nome)

        if df_OF.size == 0:
            df_OF = None
        elif (df_OF.index[0] + 1) != df_OF.index[1]: # É possível que haja uma coincidencia e repita um valor de Freq.
            #df_OF = df_OF.drop(index = df_OF.index[-1]) # Remove primeiro ponto se ele não for contínuo com o próximo
            df_OF = df_OF.drop(df_OF.index[-1])
    else:
        df_OT = None
        df_OF = None
    
    return df_CF, df_OT, df_OF

    
def main():
    nomes = [arq.split(' ')[0] for arq in glob.glob('*.txt')]
    arquivos = glob.glob('*.txt')

    for nome, arq in zip(nomes, arquivos):
        print('Tratando {0}'.format(arq))
        try:
            pd_temp = pd.read_csv(arq, delimiter=';', header=4, 
                        names=["serie", "GP", "Eta", "w", "G1", "G2", "T", "Tau", 'lixo'], encoding='latin1', decimal=',')
        except:
            pd_temp = pd.read_csv(arq, delimiter=';', header=4, 
                        names=["serie", "GP", "Tau", "GP", "Eta", "t", "t_seg", "lixo"], encoding='latin1', decimal=',')
        
        
        temp_CF, temp_OT, temp_OF = extraction(pd_temp, nome=nome + ' ')

        if type(temp_CF) != type(None):
            counter = 0
            while os.path.isfile('CF_{0}--{1}.csv'.format(nome, counter)):
                counter += 1
            temp_CF.to_csv('CF_{0}--{1}.csv'.format(nome, counter), sep=';', encoding='utf8', index=False, decimal=',')
            
        if type(temp_OT) != type(None):
            counter = 0 
            while os.path.isfile('OT_{0}--{1}.csv'.format(nome, counter)):
                counter += 1
            temp_OT.to_csv('OT_{0}--{1}.csv'.format(nome, counter), sep=';', encoding='utf8', index=False, decimal=',')
            
        if type(temp_OF) != type(None):
            counter = 0 
            while os.path.isfile('OF_{0}--{1}.csv'.format(nome, counter)):
                counter += 1
            temp_OF.to_csv('OF_{0}--{1}.csv'.format(nome, counter), sep=';', encoding='utf8', index=False, decimal=',')

        
if __name__ == '__main__':
    main()