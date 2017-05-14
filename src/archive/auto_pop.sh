for R in 75 80 85; do sed -i "s/^\(R=\).*/\1$R/" config.mk; make;done 
