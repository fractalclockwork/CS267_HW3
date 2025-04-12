export MY_DATA=$SCRATCH/my_data
export HW3_DATA=/global/cfs/cdirs/mp309/cs267-spr2025/hw3-datasets

mkdir -p $MY_DATA
lfs setstripe -c 72 -S 8M $MY_DATA 
cp -r $HW3_DATA/human-chr14-synthetic* $MY_DATA
