import argparse

parser = argparse.ArgumentParser(description='AD-MIF', formatter_class=argparse.ArgumentDefaultsHelpFormatter)

parser.add_argument('--name', type=str, default="AD")
parser.add_argument('--seed', type=int, default=0)
parser.add_argument('--rec_epoch', type=int, default=200)
parser.add_argument('--fus_epoch', type=int, default=1000)
parser.add_argument('--epoch', type=int, default=4000)
parser.add_argument('--pretrain', type=bool, default=False)

parser.add_argument('--k', type=int, default=10)
parser.add_argument('--alpha_value', type=float, default=0.1)
parser.add_argument('--lambda1', type=float, default=10)
parser.add_argument('--lambda2', type=float, default=0.1)
parser.add_argument('--lambda3', type=float, default=10)
parser.add_argument('--method', type=str, default='euc')
parser.add_argument('--first_view', type=str, default='ATAC')
parser.add_argument('--lr', type=float, default=1e-3)

parser.add_argument('--n_d1', type=int, default=4400)
parser.add_argument('--n_d2', type=int, default=4400)
parser.add_argument('--n_z', type=int, default=20)

parser.add_argument('--ae_n_enc_1', type=int, default=256)
parser.add_argument('--ae_n_enc_2', type=int, default=128)
parser.add_argument('--ae_n_dec_1', type=int, default=128)
parser.add_argument('--ae_n_dec_2', type=int, default=256)

parser.add_argument('--gae_n_enc_1', type=int, default=256)
parser.add_argument('--gae_n_enc_2', type=int, default=128)
parser.add_argument('--gae_n_dec_1', type=int, default=128)
parser.add_argument('--gae_n_dec_2', type=int, default=256)
parser.add_argument('--dropout', type=float, default=0)

parser.add_argument('--acc', type=float, default=0)
parser.add_argument('--nmi', type=float, default=0)
parser.add_argument('--ari', type=float, default=0)
parser.add_argument('--ami', type=float, default=0)

args = parser.parse_args()
