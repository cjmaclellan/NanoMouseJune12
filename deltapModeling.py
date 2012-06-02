# Read DAKOTA parameters file (aprepro or standard format) and call a
# Python module for fem analysis.

# DAKOTA will execute this script as
#   deltapModeling.py params.in results.out
# so sys.argv[1] will be the parameters file and
#    sys.argv[2] will be the results file to return to DAKOTA

# necessary python modules
import sys
import re
import os
import math

def deltapModeling(**kwargs):
  """
  treatment planning model 
  """
  # import petsc and numpy
  import petsc4py, numpy
  # init petsc
  PetscOptions =  sys.argv
  PetscOptions.append("-ksp_monitor")
  PetscOptions.append("-ksp_rtol")
  PetscOptions.append("1.0e-15")
  #PetscOptions.append("-help")
  #PetscOptions.append("-idb")
  petsc4py.init(PetscOptions)
  #
  # break processors into separate communicators
  from petsc4py import PETSc
  petscRank = PETSc.COMM_WORLD.getRank()
  petscSize = PETSc.COMM_WORLD.Get_size()
  sys.stdout.write("petsc rank %d petsc nproc %d\n" % (petscRank, petscSize))

  # set shell context
  # TODO import vtk should be called after femLibrary ???? 
  # FIXME WHY IS THIS????
  import femLibrary
  # initialize libMesh data structures
  libMeshInit = femLibrary.PyLibMeshInit(PetscOptions,PETSc.COMM_WORLD) 
  
  # store control variables
  getpot = femLibrary.PylibMeshGetPot(PetscOptions) 
  # from Duck table 2.11 pig dermis 
  getpot.SetIniValue( "material/specific_heat","3280.0" ) 
  # set ambient temperature 
  getpot.SetIniValue( "initial_condition/u_init","0.0" ) 
  # from Duck 2.7 Rat tumor measured in vivo
  getpot.SetIniValue( "thermal_conductivity/k_0_tumor","0.32" ) 
  getpot.SetIniValue( "thermal_conductivity/k_0_healthy",kwargs['cv']['k_0_healthy'] ) 
  # water properties at http://www.d-a-instruments.com/light_absorption.html
  getpot.SetIniValue( "optical/mu_a_healthy", kwargs['cv']['mu_a_healthy'] ) 
  # FIXME large mu_s (> 30) in agar causing negative fluence to satisfy BC 
  getpot.SetIniValue( "optical/mu_s_healthy",
              kwargs['cv']['mu_s_healthy'] ) 
  # from AE paper
  #http://scitation.aip.org/journals/doc/MPHYA6-ft/vol_36/iss_4/1351_1.html#F3

  #For SPIOs
  #getpot.SetIniValue( "optical/guass_radius","0.0035" ) 
  #For NR/NS
  getpot.SetIniValue( "optical/guass_radius",".003" )

  # 1-300
  getpot.SetIniValue( "optical/mu_a_tumor",
              kwargs['cv']['mu_a_tumor'] ) 
  # 1-300
  getpot.SetIniValue( "optical/mu_s_tumor",
              kwargs['cv']['mu_s_tumor'] ) 
  #getpot.SetIniValue( "optical/mu_a_tumor","71.0" ) 
  #getpot.SetIniValue( "optical/mu_s_tumor","89.0" ) 
  # .9  - .99
  getpot.SetIniValue( "optical/anfact",kwargs['cv']['anfact'] ) 
  #agar length

  #for SPIOs
  #getpot.SetIniValue( "optical/agar_length","0.0206" ) 
  #for NR/NS
  getpot.SetIniValue( "optical/agar_length","0.023" )

  getpot.SetIniValue("optical/refractive_index","1.0")
  #
  #  given the original orientation as two points along the centerline z = x2 -x1
  #     the transformed orienteation would be \hat{z} = A x2 + b - A x1 - b = A z
  #  ie transformation w/o translation which is exactly w/ vtk has implemented w/ TransformVector
  #  TransformVector = TransformPoint - the transation
  #Setup Affine Transformation for registration
  RotationMatrix = [[1.,0.,0.],
                    [0.,1.,0.],
                    [0.,0.,1.]]
  Translation =     [0.,0.,0.]
  # original coordinate system laser input
  # For SPIOs 67_11
  #laserTip         =  [0.,-.000625,.035] 
  # For SPIOs 67_10
  #laserTip         =  [0.,-.0001562,.035]
  # For SPIOs 335_11
  #laserTip         =  [0.,-.0014062,.035]
  # For NR		
  #laserTip         =  [0.,.0002,.035]
  # For NS
  laserTip         =  [0.00001,.000001,.000001]

  laserOrientation =  [0.0,0.,-1.0] 
  
  import vtk
  import vtk.util.numpy_support as vtkNumPy 
  # echo vtk version info
  print "using vtk version", vtk.vtkVersion.GetVTKVersion()
  # FIXME  notice that order of operations is IMPORTANT
  # FIXME   translation followed by rotation will give different results
  # FIXME   than rotation followed by translation
  # FIXME  Translate -> RotateZ -> RotateY -> RotateX -> Scale seems to be the order of paraview
  AffineTransform = vtk.vtkTransform()
  # should be in meters

  # For NS
  #AffineTransform.Translate([ .01320,-.0037,0.00000010])
  #AffineTransform.Translate([ .005,-.0035,0.0])
  AffineTransform.Translate([-.002,-.0022,0.0])
  AffineTransform.RotateZ( 0.0 )
  AffineTransform.RotateY( -90.0 )
  AffineTransform.RotateX( 0.0 )
  AffineTransform.Scale([1.,1.,1.])


  # get homogenius 4x4 matrix  of the form
  #               A | b
  #    matrix =   -----
  #               0 | 1
  #   
  matrix = AffineTransform.GetConcatenatedTransform(0).GetMatrix()
  #print matrix 
  RotationMatrix = [[matrix.GetElement(0,0),matrix.GetElement(0,1),matrix.GetElement(0,2)],
                    [matrix.GetElement(1,0),matrix.GetElement(1,1),matrix.GetElement(1,2)],
                    [matrix.GetElement(2,0),matrix.GetElement(2,1),matrix.GetElement(2,2)]]
  Translation =     [matrix.GetElement(0,3),matrix.GetElement(1,3),matrix.GetElement(2,3)] 
  #print RotationMatrix ,Translation 
  

  laserTip         =  AffineTransform.TransformPoint(  laserTip )
  laserOrientation =  AffineTransform.TransformVector( laserOrientation )

  # set laser orientation values
  getpot.SetIniValue( "probe/x_0","%f" % laserTip[0]) 
  getpot.SetIniValue( "probe/y_0","%f" % laserTip[1]) 
  getpot.SetIniValue( "probe/z_0","%f" % laserTip[2]) 
  getpot.SetIniValue( "probe/x_orientation","%f" % laserOrientation[0] ) 
  getpot.SetIniValue( "probe/y_orientation","%f" % laserOrientation[1] ) 
  getpot.SetIniValue( "probe/z_orientation","%f" % laserOrientation[2] ) 
  
  # initialize FEM Mesh
  femMesh = femLibrary.PylibMeshMesh()
  # must setup Ini File first
  #For SPIOs
  #femMesh.SetupUnStructuredGrid( "/work/01741/cmaclell/data/mdacc/deltap_phantom_oct10/phantomMeshFullTess.e",0,RotationMatrix, Translation  ) 
  #For NS/NR
  RotationMatrix = [[1,0,0],
                    [0,1,0],
                    [0,0,1]]
  Translation =     [.0000001,.00000001,.000001]
  #
  femMesh.SetupUnStructuredGrid( "/work/01741/cmaclell/data/mdacc/NanoMouseJune12/sphereMesh.e",0,RotationMatrix, Translation  )
  #femMesh.SetupUnStructuredGrid( "phantomMesh.e",0,RotationMatrix, Translation  )

  MeshOutputFile = "fem_data.%04d.e" % kwargs['fileID'] 
  #femMes.SetupStructuredGrid( (10,10,4) ,[0.0,1.0],[0.0,2.0],[0.0,1.0]) 
  
  # add the data structures for the background system solve
  # set deltat, number of time steps, power profile, and add system
  nsubstep = 1
  #for SPIOs
  #acquisitionTime = 4.9305
  #for NS/NR
  acquisitionTime = 5.09
  deltat = acquisitionTime / nsubstep
  ntime  = 60 
  eqnSystems =  femLibrary.PylibMeshEquationSystems(femMesh,getpot)
  #For SPIOs
  #getpot.SetIniPower(nsubstep,  [ [1,6,30,ntime],[1.0,0.0,4.5,0.0] ])
  #For NS/NR
  getpot.SetIniPower(nsubstep,  [ [1,6,42,ntime],[1.0,0.0,1.13,0.0] ])
  deltapSystem = eqnSystems.AddPennesDeltaPSystem("StateSystem",deltat) 
  deltapSystem.AddStorageVectors(ntime)

  # hold imaging
  mrtiSystem = eqnSystems.AddExplicitSystem( "MRTI" ) 
  mrtiSystem.AddFirstLagrangeVariable( "u0*" ) 
  mrtiSystem.AddStorageVectors(ntime)
  maskSystem = eqnSystems.AddExplicitSystem( "ImageMask" ) 
  maskSystem.AddFirstLagrangeVariable( "mask" ) 
  
  # initialize libMesh data structures
  eqnSystems.init( ) 
  
  # print info
  eqnSystems.PrintSelf() 
  
  # write IC
  exodusII_IO = femLibrary.PylibMeshExodusII_IO(femMesh)
  exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, 1, 0.0 )  
  
  # read imaging data geometry that will be used to project FEM data onto
  #For 67_11
  #vtkReader.SetFileName('/work/01741/cmaclell/data/mdacc/deltap_phantom_oct10/spio/spioVTK/67_11/tmap_67_11.0000.vtk')
  #For 67_10
  #vtkReader.SetFileName('/work/01741/cmaclell/data/mdacc/deltap_phantom_oct10/spio/spioVTK/67_10/tmap_67_10.0000.vtk')
  #For 335_11 
  #vtkReader.SetFileName('/work/01741/cmaclell/data/mdacc/deltap_phantom_oct10/spio/spioVTK/335_11/tmap_335_11.0000.vtk')
  #For NS
  #vtkReader.SetFileName('/work/01741/cmaclell/data/mdacc/deltap_phantom_oct10/nrtmapsVTK/S695/S695.0000.vtk') 
  #For NR
  imageFileNameTemplate = '/work/01741/cmaclell/data/mdacc/NanoMouseJune12/matlab_VTK/control_1_tmap.%04d.vtk'
  #imageFileNameTemplate = "/share/work/fuentes/deltap_phantom_oct10/nrtmapsVTK/R695/R695.%04d.vtk"
  #imageFileNameTemplate = "/data/fuentes/mdacc/deltap_phantom_oct10/nrtmapsVTK/R695/R695.%04d.vtk"

  #vtkReader = vtk.vtkXMLImageDataReader() 
  vtkReader = vtk.vtkDataSetReader() 
  vtkReader.SetFileName( imageFileNameTemplate % 0 )
  vtkReader.Update()
  templateImage = vtkReader.GetOutput()
  dimensions = templateImage.GetDimensions()
  spacing = templateImage.GetSpacing()
  origin  = templateImage.GetOrigin()
  print spacing, origin, dimensions
  femImaging = femLibrary.PytttkImaging(getpot, dimensions ,origin,spacing) 
  
  ObjectiveFunction = 0.0
  # loop over time steps and solve
  for timeID in range(1,ntime*nsubstep):
  #for timeID in range(1,10):
     # project imaging onto fem mesh
     vtkImageReader = vtk.vtkDataSetReader() 
     vtkImageReader.SetFileName( imageFileNameTemplate % timeID )
     vtkImageReader.Update() 
     image_cells = vtkImageReader.GetOutput().GetPointData() 
     data_array = vtkNumPy.vtk_to_numpy(image_cells.GetArray('scalars')) 
     v1 = PETSc.Vec().createWithArray(data_array, comm=PETSc.COMM_SELF)
     femImaging.ProjectImagingToFEMMesh("MRTI",0.0,v1,eqnSystems)  
     mrtiSystem.StoreSystemTimeStep(timeID ) 
  
     # create image mask 
     #  The = operator for numpy arrays just copies by reference
     #   .copy() provides a deep copy (ie physical memory copy)
     image_mask = data_array.copy().reshape(dimensions,order='F')
     # Set all image pixels to large value
     largeValue = 1.e6
     image_mask[:,:] = 1 
     # RMS error will be computed within this ROI/VOI imagemask[xcoords/column,ycoords/row]
     #image_mask[93:153,52:112] = 1.0
     #image_mask[98:158,46:106] = 1.0
     v2 = PETSc.Vec().createWithArray(image_mask, comm=PETSc.COMM_SELF)
     femImaging.ProjectImagingToFEMMesh("ImageMask",largeValue,v2,eqnSystems)  
     #print mrti_array
     #print type(mrti_array)

     print "time step = " ,timeID
     eqnSystems.UpdatePetscFEMSystemTimeStep("StateSystem",timeID ) 
     deltapSystem.SystemSolve()
     #eqnSystems.StoreTransientSystemTimeStep("StateSystem",timeID ) 
  
     # accumulate objective function
     # 
     # qoi  = ( (FEM - MRTI) / ImageMask)^2
     # 
     qoi = femLibrary.WeightedL2Norm( deltapSystem,"u0",
                                      mrtiSystem,"u0*",
                                      maskSystem,"mask" ) 
     ObjectiveFunction =( ObjectiveFunction + qoi) 
    
     # control write output
     writeControl = True
     if ( timeID%nsubstep == 0 and writeControl ):
       # write exodus file
       exodusII_IO.WriteTimeStep(MeshOutputFile,eqnSystems, timeID+1, timeID*deltat )  
       # Interpolate FEM onto imaging data structures
       if (vtk != None):
         vtkExodusIIReader = vtk.vtkExodusIIReader()
         vtkExodusIIReader.SetFileName(MeshOutputFile )
         vtkExodusIIReader.SetPointResultArrayStatus("u0",1)
         vtkExodusIIReader.SetTimeStep(timeID-1) 
         vtkExodusIIReader.Update()
         # reflect
         vtkReflect = vtk.vtkReflectionFilter()
         vtkReflect.SetPlaneToYMax()
         vtkReflect.SetInput( vtkExodusIIReader.GetOutput() )
         vtkReflect.Update()
         # reuse ShiftScale Geometry
         vtkResample = vtk.vtkCompositeDataProbeFilter()
         vtkResample.SetInput( vtkImageReader.GetOutput() )
         vtkResample.SetSource( vtkReflect.GetOutput() ) 
         vtkResample.Update()
         fem_point_data= vtkResample.GetOutput().GetPointData() 
         fem_array = vtkNumPy.vtk_to_numpy(fem_point_data.GetArray('u0')) 
         #print fem_array 
         #print type(fem_array )
       # only rank 0 should write
       #if ( petscRank == 0 ):
       #   print "writing ", timeID
       #   vtkTemperatureWriter = vtk.vtkDataSetWriter()
       #   vtkTemperatureWriter.SetFileTypeToBinary()
       #   vtkTemperatureWriter.SetFileName("invspio_67_11.%04d.vtk" % timeID )
       #   vtkTemperatureWriter.SetInput(vtkResample.GetOutput())
       #   vtkTemperatureWriter.Update()
  print 'Objective Fn'
  print ObjectiveFunction
  retval = dict([])
  retval['fns'] = [ObjectiveFunction *10000000.]
  retval['rank'] = petscRank 
  return(retval)
# end def deltapModeling(**kwargs):
########################################################################################

# ----------------------------
# Parse DAKOTA parameters file
# ----------------------------

# setup regular expressions for parameter/label matching
e = '-?(?:\\d+\\.?\\d*|\\.\\d+)[eEdD](?:\\+|-)?\\d+' # exponential notation
f = '-?\\d+\\.\\d*|-?\\.\\d+'                        # floating point
i = '-?\\d+'                                         # integer
value = e+'|'+f+'|'+i                                # numeric field
tag = '\\w+(?::\\w+)*'                               # text tag field

# regular expression for aprepro parameters format
aprepro_regex = re.compile('^\s*\{\s*(' + tag + ')\s*=\s*(' + value +')\s*\}$')
# regular expression for standard parameters format
standard_regex = re.compile('^\s*(' + value +')\s+(' + tag + ')$')

# open DAKOTA parameters file for reading
paramsfile = open(sys.argv[1], 'r')
fileID = int(sys.argv[1].split(".").pop())

# extract the parameters from the file and store in a dictionary
paramsdict = {}
for line in paramsfile:
    m = aprepro_regex.match(line)
    if m:
        paramsdict[m.group(1)] = m.group(2)
    else:
        m = standard_regex.match(line)
        if m:
            paramsdict[m.group(2)] = m.group(1)

paramsfile.close()

# crude error checking; handle both standard and aprepro cases
num_vars = 0
if ('variables' in paramsdict):
    num_vars = int(paramsdict['variables'])
elif ('DAKOTA_VARS' in paramsdict):
    num_vars = int(paramsdict['DAKOTA_VARS'])

num_fns = 0
if ('functions' in paramsdict):
    num_fns = int(paramsdict['functions'])
elif ('DAKOTA_FNS' in paramsdict):
    num_fns = int(paramsdict['DAKOTA_FNS'])

# -------------------------------
# Convert and send to application
# -------------------------------

# set up the data structures the rosenbrock analysis code expects
# for this simple example, put all the variables into a single hardwired array
continuous_vars = { 
                    'k_0_healthy' :paramsdict['k_0_healthy'  ],
                    'k_0_tumor'   :'.63' ,
                    'mu_a_healthy':paramsdict['mu_a_healthy'  ],
                    'mu_a_tumor'  :'1.0'   
                  }

try:
   continuous_vars['w_0_healthy'] = paramsdict['w_0_healthy' ]  
   continuous_vars['w_0_tumor'  ] = paramsdict['w_0_tumor'   ] 
except KeyError:
   continuous_vars['w_0_healthy'] = "0.0"
   continuous_vars['w_0_tumor'  ] = "0.0"

try:
   continuous_vars['anfact'] = paramsdict['anfact'   ] 
except KeyError:
   continuous_vars['anfact'] = "0.9"

try:
   continuous_vars['mu_s_healthy'] = paramsdict['mu_s_healthy']
   continuous_vars['mu_s_tumor'  ] = paramsdict['mu_s_tumor'  ]
except KeyError:
   #anfact       = float(continuous_vars['anfact'] )
   #od_healthy   = float('.105')
   #od_tumor     = float('.707')
   mu_s = float('310')
   #mu_a_tumor   = float('2')
   #Mutr=ln(10)*OD/.01  #  .01 --> in meters  
   #mu_s = (mutr-mua)/(1-g)
   #mu_tr_healthy= math.log(10) * od_healthy / 0.01
   #mu_tr_tumor  = math.log(10) * od_tumor   / 0.01
   continuous_vars['mu_s_healthy'] = "%f" % mu_s
   continuous_vars['mu_s_tumor'  ] = "%f" % mu_s

try:
   continuous_vars['x_translate'] = float( '0.0000001' )
except KeyError:
   continuous_vars['x_translate'] = -0.0055

try:
  active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
except KeyError:
  active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 

# set a dictionary for passing to rosenbrock via Python kwargs
fem_params              = {}
fem_params['cv']        = continuous_vars
fem_params['asv']       = active_set_vector
fem_params['functions'] = num_fns
fem_params['fileID']    = fileID 

# execute the rosenbrock analysis as a separate Python module
print "Running deltap model..."
fem_results = deltapModeling(**fem_params)
print "deltap complete."
print continuous_vars


# ----------------------------
# Return the results to DAKOTA
# ----------------------------

if (fem_results['rank'] == 0 ):
  # write the results.out file for return to DAKOTA
  # this example only has a single function, so make some assumptions;
  # not processing DVV
  outfile = open('results.out.tmp.%d' % fileID, 'w')
  
  # write functions
  for func_ind in range(0, num_fns):
      if (active_set_vector[func_ind] & 1):
          functions = fem_results['fns']    
          outfile.write(str(functions[func_ind]) + ' f' + str(func_ind) + '\n')
  
  ## write gradients
  #for func_ind in range(0, num_fns):
  #    if (active_set_vector[func_ind] & 2):
  #        grad = rosen_results['fnGrads'][func_ind]
  #        outfile.write('[ ')
  #        for deriv in grad: 
  #            outfile.write(str(deriv) + ' ')
  #        outfile.write(']\n')
  #
  ## write Hessians
  #for func_ind in range(0, num_fns):
  #    if (active_set_vector[func_ind] & 4):
  #        hessian = rosen_results['fnHessians'][func_ind]
  #        outfile.write('[[ ')
  #        for hessrow in hessian:
  #            for hesscol in hessrow:
  #                outfile.write(str(hesscol) + ' ')
  #            outfile.write('\n')
  #        outfile.write(']]')
  #
  outfile.close();outfile.flush
  #
  # move the temporary results file to the one DAKOTA expects
  import shutil
  shutil.move('results.out.tmp.%d' % fileID, sys.argv[2])
  #os.system('mv results.out.tmp ' + sys.argv[2])
