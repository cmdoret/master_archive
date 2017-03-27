# Just a small script to rename samples by adding their sex.
import pandas as pd

# Loading data files
clutch = pd.read_csv("F4_table.csv")
adap = pd.read_csv("adaptors_radwasp7.csv")

sex = {}
for i,j in clutch.iterrows():
    sex[j['Number']] = j['Sex']

for n in range(adap.shape[0]):
    try:
        adap.loc[(n,'sample')] += ('_' + sex[adap.loc[(n,'sample')]])
    except KeyError:
        print("contaminated sample")
    print(adap.loc[(n,'sample')])

adap.to_csv("new_adaptors.csv",index=False)
