import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler
import xgboost as xgb
from sklearn.ensemble import RandomForestRegressor

def prepare_data():
    import numpy as np



    s_data = np.load('snp.npy')
    s_data[np.isnan(s_data)] = -1
    labels = np.load('pheno.npy')

    # Concatenate the labels and s_data horizontally
    combined_data = np.hstack((labels, s_data))

    # Write the data to a CSV file
    with open('data_bridge.csv', 'w') as file:
        # Write the header
        file.write("ID,label,snp\n")

        for i, row in enumerate(combined_data, start=1):
            # Convert nan to 'np.nan' and other values to their string representation
            values = ["np.nan" if np.isnan(val) else str(val) for val in row]

            # Join the values with a comma and write to the file
            file.write("{},{}\n".format(i, '\t'.join(values)))
def run_tree():
    snp_data = np.load('snp.npy')
    print (snp_data)
    snp_data[np.isnan(snp_data)] = -1
    phenotype_labels = np.load('pheno.npy')

    #print (phenotype_labels)
    # Split data into training and testing sets
    X_train, X_test, y_train, y_test = train_test_split(snp_data, phenotype_labels, test_size=0.2, random_state=42)

    # Standardize the input features (optional but recommended)
    scaler = StandardScaler()
    X_train = scaler.fit_transform(X_train)
    X_test = scaler.transform(X_test)

    # Define and train an XGBoost model
    # xgboost1
    #xgb_model = xgb.XGBRegressor(n_estimators=100, learning_rate=0.1, max_depth=3, random_state=42)

    #xgb_model = RandomForestRegressor(n_estimators=100, random_state=42)
    # xgboost2
    xgb_model = xgb.XGBRegressor(n_estimators=25
                                      ,
                                      max_depth=20,
                                      min_child_weight=10,
                                      subsample=0.8,
                                      eta=0.3,
                                      colsample_bytree=0.7)
    # Fit the model to the training data
    xgb_model.fit(X_train, y_train)

    # Evaluate the model on the test set
    test_predictions = xgb_model.predict(X_test)
    test_loss = np.mean((test_predictions - y_test)**2)  # Mean squared error
    print(f'Test Loss: {test_loss}')
    f = open('out_xgboos_bridge2.txt', 'w')

    #feature_list = open('C:/gwas_test_data/test/' + "WGS300_005_0020.bim", 'r')
    #feature_list = [line.rstrip() for line in feature_list.readlines()]
    df = pd.read_csv('C:/gwas_test_data/test/' + "WGS300_005_0020.bim", delimiter='\t')
    feature_list = df.iloc[:, 1].tolist()
    print (feature_list)
    for col, score in zip(feature_list, xgb_model.feature_importances_):
        f.write(str(col) + ' ' + str(score) + '\n')


def plot():
    #df = pd.read_csv('out_xgboos_bridge2.txt', delimiter=' ', header=None)
    df = pd.read_csv('C:/Users/lueck/PycharmProjects/genbank2/test/_out_nn_bridge100.txt', delimiter=' ', header=None)

    df.columns = ['snp', 'value']
    #print (df.dtypes)
    #print(df)
    #df = df.sort_values(by=['value'], ascending=False)
    #ax1 = df.plot.scatter(x='value',
                        #  y='snp')
    #plt.savefig('out_xgboos_bridge.png')
    feature_list = df['value'].to_numpy()
    plt.plot(feature_list)
    plt.show()
    #print (feature_list)
    #plt.axhline(y=feature_list, color='b', linestyle='--')
    # plt.ylabel('saliency value', fontdict=None, labelpad=None, fontsize=15)
    #
    # plt.xlabel('SNPs', fontdict=None, labelpad=None, fontsize=15)
    #
    #plt.savefig('sal.png')
    #

    # print (len(xgb_model.feature_importances_))
    # xgb.plot_importance(xgb_model, max_num_features=20)
    # plt.rcParams['figure.figsize'] = [50, 50]
    # plt.savefig('rf.png')

#run_tree()
plot()
