USAGE <MODEL_PATH> <OUTPUT_PATH> <TARGET_ID> <DNCON Feature path> <INTRACHAIN Feature path> <SS8 Feature path> <trRoseetta Feature path>
e.g
python DRCON_pred.py /gpfs/alpine/proj-shared/bif132/raj/codes/pytroch_codes/updated_history/weighths_82 /gpfs/alpine/proj-shared/bif132/raj/codes/pytroch_codes/ 3GWR /gpfs/alpine/proj-shared/bif132/raj/dimer/Deephomo_data/dncon_feat_only/feat-3GWR.txt /gpfs/alpine/proj-shared/bif132/raj/dimer/Deephomo_data/intra_cmap/3GWR.cmap /gpfs/alpine/proj-shared/bif132/raj/dimer/Deephomo_data/ss8_one_hot/3GWR.feat_ss8 /gpfs/alpine/proj-shared/bif132/raj/dimer/Deephomo_data/tr_features/3GWR.npz


##TESTED ON THE FOLLOWING LIBRARY PACKAGE##
python                    3.8.8                h836d2c2_4
pytorch                   1.7.1           cuda10.2_py38_2
pytorch-base              1.7.1           cuda10.2_py38_8
torchtext                 0.8.1                    py38_4
torchvision               0.8.2           cuda10.2_py38_0
torchvision-base          0.8.2           cuda10.2_py38_6
