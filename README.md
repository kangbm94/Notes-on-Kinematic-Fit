# Notes-on-Kinematic-Fit -- provided with example codes.
# author: Byungmin Kang, email: kangbmw2@naver.com 
#      -----------
#      | KinFit  |  <- Abstract class for KF procedure. User does not need modify.
#      -----------
#           |
#      ---------------
#      V             V
#   ----------   ----------
#   |Cas.Fit |   |F.V.Fit |   <- Concrete class for application. User should define the details, depending on their situations.
#   ----------   ----------
# 
#  User can inherit KinFit class to define their own class for application. Variables of interest, their constraints, derivatives, and covariance matrix should be defined by hand.
#  Process related to optimization is dealt in KinFit class. Once user defined their applications properly, statistical variables like chi2, p-values, or pull distrubutions can be extracted from KinFit.


