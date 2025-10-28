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
# Follwoings are AI Generated description #
+-------------------------------------------------------------+
|                    KinematicFitter Core                     |
+-------------------------------------------------------------+
|Setup phase                                                  |
|   SetVariance(var)      [src/KinFit.cc:15]                  |
|       └── seeds Variancies/ScalingMats                      |
|   AddOffdiagonals(Cov)  [src/KinFit.cc:58]                  |
|       └── adjusts Variancies[0] to stay positive definite   |
|   RotateVariance(J)     [src/KinFit.cc:461] (optional)      |
+-------------------------------------------------------------+
|Iterative solve                                              |
|   DoKinematicFit(Do)    [src/KinFit.cc:393]                 |
|       ├── guard: Finalize() if Do=false                     |
|       ├── loop: ProcessStep() until χ² stabilizes/MAX step  |
|       └── Finalize() to pick best step & store pulls        |
|         (Finalize uses FMats/Pulls to choose best solution) |
+-------------------------------------------------------------+
|Single step mechanics                                        |
|   ProcessStep()          [src/KinFit.cc:109]                |
|       ├── load current Measurements/Unknowns                |
|       ├── call virtual SetConstraints() (user hook)         |
|       ├── build residual r, covariance s, multipliers λ     |
|       ├── update Unknowns/Measurements/Variancies           |
|       ├── collect χ², pulls, covariances for this step      |
|       └── call virtual SampleStepPoint(step) (user hook)    |
+-------------------------------------------------------------+
|Housekeeping                                                 |
|   Clear()               [src/KinFit.cc:425] resets state    |
|   TransposeMatrix(M)    [src/KinFit.cc:448] helper          |
+-------------------------------------------------------------+
