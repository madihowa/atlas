from CSV_read_write import *
import os
dir = "Results_2020-07-07_run_2"

os.chdir(dir)
from DNN import Make_network, compile_NN
os.chdir("../")



def AddEM():


    test, train = read_csv_file()

    list_input = ReadInputVaribles()
    file = find_min_weights(dir)
    NN_model = Make_network(test[list_input].shape[1])

    NN_model.load_weights(file)

    NN_model = compile_NN(NN_model)

    test['EM_Char'] = NN_model.predict(test[list_input])
    train['EM_Char'] = NN_model.predict(train[list_input])

    write_csv_file(test, "piEMpro", "test.csv")
    write_csv_file(train, "piEMpro", "train.csv")


AddEM()
