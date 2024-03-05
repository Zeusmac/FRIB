import numpy as np

## p2*exp(p1*x) + bkgrond

def half_life(p1,error):
    
    return [np.log(2)/np.abs(p1) *1e-6,np.abs(np.log(2)/pow(p1,2))*error*1e-6]

p1Na31 = float(6.62084e-08)
bkground = float(3.25818e+01)
bkgrounderror = float(3.12795e-01) 
p1errorNa31 = float(8.52287e-09)
p2errorNa31 = float(4.62185e+00)

tNa31 = np.log(2)/np.abs(p1Na31)
errorNa31 = np.abs(np.log(2)/pow(p1Na31,2))*p1errorNa31

print("Na31 half life = %5f" % (tNa31*1.0e-6))
print("Na31 error = %5f" % (errorNa31*1.0e-6))## p2*exp(p1*x) + bkgrond

##p2Na32 = float(2.54627e+01)
p2Na32 = float(4.66075e-08)
##p1Na32 = float(-4.66073e-08)

print("Na32 Halflife,error:", half_life(1.11273e-07,2.17910e-08 ))
print("Daughter Halflife,error:", half_life(2.03737e-03 ,1.85181e-02))
print(half_life(2.27018e-05,7.36964e-05))

p1errorNa32 = float(5.54352e-09)
p2errorNa32 = float(4.62185e+00)

# tNa32 = np.log(2)/np.abs(p1Na32)
# errorNa32 = np.abs(np.log(2)/pow(p1Na32,2))*p1errorNa32

# print("Na32 half life = %5f" % (tNa32*1.0e-6))
# print("Na32 error = %5f" % (errorNa32*1.0e-6))

p1Na33 = float(-1.38650e-07 )

p1errorNa33 = float(2.52719e-08)


tNa33 = np.log(2)/np.abs(p1Na33)
errorNa33 = np.abs(np.log(2)/pow(p1Na33,2))*p1errorNa33

print("Na33 half life = %5f" % (tNa33*1.0e-6))
print("Na33 error = %5f" % (errorNa33*1.0e-6))

p1Na30 = float(-6.82138e-08)

p1errorNa30 = float(2.00402e-08)


tNa30 = np.log(2)/np.abs(p1Na30)
errorNa30 = np.abs(np.log(2)/pow(p1Na30,2))*p1errorNa30

print("Na30 half life = %5f" % (tNa30*1.0e-6))
print("Na30 error = %5f" % (errorNa30*1.0e-6))

print("Na31 Halflife,error:", half_life(6.62084e-08,8.52287e-09))
print("daughter Na31 Halflife,error:", half_life(2.95564e-03,7.58961e-03))








