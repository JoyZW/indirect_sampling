# Regarding the unbiasedness of the sample-dependent strategy, let me make the same suggestion again. 
# 
# - Take the small example in my note, and assume SRS of 2 customers initially. 
# There are 6 (4 choose 2) distinct initial (customer) samples, all with the same probability 1/6.
# 
# - Fix any sequence of the 6 samples. 
# The true expectation of the estimator (2) under the sample-dependent strategy, 
# given that sequence, can be obtained numerically by going over all these 6 samples 
# one after another (which is equivalent to repeated sampling infinitely many times, since p(s)= 1/6).
# 
# i.e. E(muhatbar) = 1/6*estimate on sample 1 + .... + 1/6*estimate on sample 6

# - Now, actually go through the 6 initial samples one after another. 
# * given each initial sample, determine the additional edges to H_|3^S according
# to my description, if these exist; (NB. The whole edge set H_|3^S will most 
# likely be fixed before you have gone through all the 6 samples.)
# * calculate the estimator (2) according to H_|3^S, which would only require the 
# part that has been fixed so far. (NB. It may be helpful to write down 
# explicitly all the relevant intermediate results at each step for yourself.)
# 
# In this way, you should be able to verify the unbiasedness for yourself. 
# I imagine that afterwards you would be able to supply a formal argument for unbiasedness.