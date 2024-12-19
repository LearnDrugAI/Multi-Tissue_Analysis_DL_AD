import numpy as np
import tqdm
from torch.optim import Adam
from time import time

from utils import *
from encoder import *
from AD_MIF import AD_MIF
from data_loader import load_data

def train(model, X_train1, A_train1, X_train2, A_train2, y_train,
          X_val1, A_val1, X_val2, A_val2, y_val, X_test1, A_test1, X_test2, A_test2, y_test,Xr,Xa,Ar,Aa,train_idx,val_idx,test_idx):

    print("Training...")

    optimizer = Adam(model.parameters(), lr=(opt.args.lr))
    # y_train = y_train-1
    # y_val = y_val-1
    # y_test = y_test-1

    y_train = torch.tensor(y_train, dtype=torch.long).to('cuda:0')
    y_val = torch.tensor(y_val, dtype=torch.long).to('cuda:0')
    y_test = torch.tensor(y_test, dtype=torch.long).to('cuda:0')

    best_auc = 0
    best_acc = 0
    best_epoch = 0
    best_model_state = None

    pbar = tqdm.tqdm(range(opt.args.epoch), ncols=200)
    for epoch in pbar:
        model.train()
        X_hat1, Z_hat1, A_hat1, X_hat2, Z_hat2, A_hat2, Q1, Q2, Z1, Z2, cons, output = model(Xr, Ar, Xa, Aa)

        criterion = nn.CrossEntropyLoss()
        loss_output = criterion(output[train_idx], y_train)

        L_DRR = drr_loss(cons)

        L_REC1 = reconstruction_loss(Xr, Ar, X_hat1, Z_hat1, A_hat1)
        L_REC2 = reconstruction_loss(Xa, Aa, X_hat2, Z_hat2, A_hat2)


        L_KL1 = distribution_loss(Q1, target_distribution(Q1[0].data))
        L_KL2 = distribution_loss(Q2, target_distribution(Q1[0].data))

        loss = L_REC1 +L_REC2+ L_DRR + opt.args.lambda3 * (L_KL1 + L_KL2) + loss_output


        optimizer.zero_grad()
        loss.backward()
        optimizer.step()

        _, predicted = torch.max(output[train_idx], 1)
        train_acc = (predicted == y_train).sum().item() / y_train.size(0)

        model.eval()
        with torch.no_grad():
            _, _, _, _, _, _, _, _, _, _, _, output = model(Xr, Ar, Xa, Aa)

            output = output[val_idx]
            val_loss_output = criterion(output, y_val)

            _, predicted_val = torch.max(output, 1)
            val_acc = (predicted_val == y_val).sum().item() / y_val.size(0)

            all_labels = y_val.cpu().numpy()
            all_predictions = torch.softmax(output, dim=1).cpu().numpy()
            all_predictions = np.array(all_predictions)[:,1]

            from sklearn.metrics import roc_auc_score


            auc_val = roc_auc_score(all_labels, all_predictions)

        if auc_val > best_auc:
            best_auc = auc_val
            best_acc = val_acc
            best_epoch = epoch
            best_model_state = model.state_dict()

        pbar.set_postfix({'train_loss': '{0:1.4f}'.format(loss.item()),
                          'train_ACC': '{0:1.4f}'.format(train_acc),
                          'val_ACC': '{0:1.4f}'.format(val_acc),
                          'val_AUC': '{0:1.4f}'.format(auc_val)})

    pbar.close()

    print(f"Best_epoch: {best_epoch}, ACC: {best_acc:.4f}, AUC: {best_auc:.4f}")

    model.load_state_dict(best_model_state)
    model.eval()
    with torch.no_grad():
        _, _, _, _, _, _, _, _, _, _, _, output = model(Xr, Ar, Xa, Aa)

        _, predicted_test = torch.max(output[test_idx], 1)
        test_acc = (predicted_test == y_test).sum().item() / y_test.size(0)

        all_labels_test = y_test.cpu().numpy()
        all_predictions_test = torch.softmax(output[test_idx], dim=1).cpu().numpy()

        all_predictions_test = np.array(all_predictions_test)[:,1]

        # auc_test = roc_auc_score(all_labels_test, all_predictions_test, multi_class='ovr')
        auc_test = roc_auc_score(all_labels_test, all_predictions_test)


    print(f"Test Results -> ACC: {test_acc:.4f}, AUC: {auc_test:.4f}")

    np.save('seed{}_label_test.npy'.format(opt.args.name, opt.args.seed), predicted_test.cpu().numpy())
    np.save('seed{}_output_test.npy'.format(opt.args.name, opt.args.seed), output[test_idx].cpu().detach().numpy())


if __name__ == '__main__':
    print("setting:")

    setup_seed(opt.args.seed)

    opt.args.device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    
    print("------------------------------")
    print("dataset       : {}".format(opt.args.name))
    print("device        : {}".format(opt.args.device))
    print("random seed   : {}".format(opt.args.seed))
    print("lambda1 value : {}".format(opt.args.lambda1))
    print("lambda2 value : {}".format(opt.args.lambda2))
    print("lambda3 value : {}".format(opt.args.lambda3))
    print("alpha value   : {:.0e}".format(opt.args.alpha_value))
    print("k value       : {}".format(opt.args.k))
    print("learning rate : {:.0e}".format(opt.args.lr))
    print("------------------------------")

    Xr, y, Ar = load_data(opt.args.name, 'gen_exp_low_dim', opt.args.method, opt.args.k,'labels_apoe', show_details=False)

    Xa, y, Aa = load_data(opt.args.name, 'gen_phy_low_dim', opt.args.method, opt.args.k,'labels_apoe', show_details=False)

    opt.args.n_clusters = int(max(y) - min(y) + 1)

    from sklearn.model_selection import train_test_split

    indices = np.arange(Xr.shape[0])

    Xr_train, Xr_temp, Ar_train, Ar_temp, Xa_train, Xa_temp, Aa_train, Aa_temp, y_train, y_temp, train_idx, temp_idx = train_test_split(
        Xr, Ar, Xa, Aa, y, indices, test_size=0.3, random_state=42,
    )

    Xr_val, Xr_test, Ar_val, Ar_test, Xa_val, Xa_test, Aa_val, Aa_test, y_val, y_test, val_idx, test_idx = train_test_split(
        Xr_temp, Ar_temp, Xa_temp, Aa_temp, y_temp, temp_idx, test_size=0.33, random_state=42,

    )

    Xr_train = numpy_to_torch(Xr_train).to(opt.args.device)
    Xr_val = numpy_to_torch(Xr_val).to(opt.args.device)
    Xr_test = numpy_to_torch(Xr_test).to(opt.args.device)

    Ar_train = numpy_to_torch(Ar_train, sparse=True).to(opt.args.device)
    Ar_val = numpy_to_torch(Ar_val, sparse=True).to(opt.args.device)
    Ar_test = numpy_to_torch(Ar_test, sparse=True).to(opt.args.device)

    Xa_train = numpy_to_torch(Xa_train).to(opt.args.device)
    Xa_val = numpy_to_torch(Xa_val).to(opt.args.device)
    Xa_test = numpy_to_torch(Xa_test).to(opt.args.device)

    Ar_train = numpy_to_torch(Aa_train, sparse=True).to(opt.args.device)
    Aa_val = numpy_to_torch(Aa_val, sparse=True).to(opt.args.device)
    Aa_test = numpy_to_torch(Aa_test, sparse=True).to(opt.args.device)

    Xr = numpy_to_torch(Xr).to(opt.args.device)
    Ar = numpy_to_torch(Ar, sparse=True).to(opt.args.device)

    Xa = numpy_to_torch(Xa).to(opt.args.device)
    Aa = numpy_to_torch(Aa, sparse=True).to(opt.args.device)

    train_idx = torch.LongTensor(train_idx).to(opt.args.device)
    val_idx = torch.LongTensor(val_idx).to(opt.args.device)
    test_idx = torch.LongTensor(test_idx).to(opt.args.device)


    ae1 = AE(
        ae_n_enc_1=opt.args.ae_n_enc_1, ae_n_enc_2=opt.args.ae_n_enc_2,
        ae_n_dec_1=opt.args.ae_n_dec_1, ae_n_dec_2=opt.args.ae_n_dec_2,
        n_input=opt.args.n_d1, n_z=opt.args.n_z).to(opt.args.device)

    ae2 = AE(
        ae_n_enc_1=opt.args.ae_n_enc_1, ae_n_enc_2=opt.args.ae_n_enc_2,
        ae_n_dec_1=opt.args.ae_n_dec_1, ae_n_dec_2=opt.args.ae_n_dec_2,
        n_input=opt.args.n_d2, n_z=opt.args.n_z).to(opt.args.device)
    
    if opt.args.pretrain:
        opt.args.dropout = 0.4
    gae1 = IGAE(
        gae_n_enc_1=opt.args.gae_n_enc_1, gae_n_enc_2=opt.args.gae_n_enc_2,
        gae_n_dec_1=opt.args.gae_n_dec_1, gae_n_dec_2=opt.args.gae_n_dec_2,
        n_input=opt.args.n_d1, n_z=opt.args.n_z, dropout=opt.args.dropout).to(opt.args.device)

    gae2 = IGAE(
        gae_n_enc_1=opt.args.gae_n_enc_1, gae_n_enc_2=opt.args.gae_n_enc_2,
        gae_n_dec_1=opt.args.gae_n_dec_1, gae_n_dec_2=opt.args.gae_n_dec_2,
        n_input=opt.args.n_d2, n_z=opt.args.n_z, dropout=opt.args.dropout).to(opt.args.device)


    t0 = time()
    model = AD_MIF(ae1, ae2, gae1, gae2, n_node=Xr.shape[0]).to(opt.args.device)

    train(model, Xr_train, Ar_train, Xa_train, Aa_train, y_train,
          Xr_val, Ar_val, Xa_val, Aa_val, y_val, Xr_test, Ar_test, Xa_test, Aa_test, y_test,Xr,Xa,Ar,Aa,train_idx,val_idx,test_idx)


    t1 = time()
    print("Time_cost: {}".format(t1-t0))
