# clarox -nographics -nojournal -noecho -batch journalfilename.py
# from cubit cmd script editor
#   f = file("file.py"); cmdtext = f.read() ; exec(cmdtext)

def SphereMesh(sphereRadius):
  """
  Vessel Radius and Distance input in mm
  """
  cubit.cmd('set developer on')
  cubit.cmd('   reset')
  
  #create sphere octant
  cubit.cmd('sphere Radius %f xpos ypos zpos' % sphereRadius)
  cubit.cmd('volume 1 scale 1 1 1')
  cubit.cmd('volume 1 name "mouseleg" ')
  idhealthy = cubit.get_id_from_name('mouseleg')
  print "id" ,idhealthy
  cubit.cmd('volume 1   copy reflect x ')
  cubit.cmd('volume all copy reflect y ')
  # imprint and merge
  cubit.cmd('imprint volume all ')
  cubit.cmd('merge   volume all ')
  
  # set size
  nelem = 6
  meshSize = sphereRadius / nelem 
  cubit.cmd('volume  all size %f' %   meshSize )
  ## # mesh 
  cubit.cmd('list volume  1')
  cubit.cmd('volume  all scheme tetprimitive')
  cubit.cmd('mesh volume  all')
  
  # export in pieces
  cubit.cmd('reset genesis')
  cubit.cmd('block 1 volume  1 2 3 4')
  cubit.cmd('block 1 name "mouseleg"  ')
  # add BC
  #cubit.cmd('skin volume all make sideset 2')
  cubit.cmd('sideset 2 surface 2 6 10 14')
  cubit.cmd('sideset 2 name "neumann" ')
  cubit.cmd('sideset 4 surface 4 8 12 16')
  cubit.cmd('sideset 4 name "fluence" ')
  #cubit.cmd('nodeset 1 volume 17 18 19')
  #cubit.cmd('nodeset 1 name "dirichletApplicator"')
  #cubit.cmd('nodeset 2 volume  9 15 20')
  #cubit.cmd('nodeset 2 name "dirichletVessel"')
  #
  # scale from [mm] to [m] and write'
  cubit.cmd('volume all scale 0.001')
  #rotate mesh to align with imaging data
  #cubit.cmd('volume all rotate 0 about x')
  #cubit.cmd('volume all rotate -90 about y')
  #cubit.cmd('volume all rotate -45 about z')
  #translate mesh to align with imaging data
  #cubit.cmd('volume all move X .01320 Y -.0037 Z 0')
  cubit.cmd('export mesh "sphereMesh.e" overwrite' )
# end def SphereMesh
##################################################################
def ParseDakotaFile(param_file):
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
  paramsfile = open(param_file, 'r')
  
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
                      'vessel_distance' :paramsdict['vessel_distance'],
                      'vessel_diameter' :paramsdict['vessel_diameter'],
                    }
  
  try:
    active_set_vector = [ int(paramsdict['ASV_%d:response_fn_%d' % (i,i) ]) for i in range(1,num_fns+1)  ] 
  except KeyError:
    active_set_vector = [ int(paramsdict['ASV_%d:obj_fn' % (i) ]) for i in range(1,num_fns+1)  ] 
  
  # set a dictionary for passing to rosenbrock via Python kwargs
  fem_params              = {}
  fem_params['cv']        = continuous_vars
  fem_params['asv']       = active_set_vector
  fem_params['functions'] = num_fns

  return fem_params              
# end def ParseDakotaFile:
##################################################################
#params = ParseDakotaFile("/data/fuentes/utsa/vasculature_july10/vessel/realization.1/pce.in")
#vessel_distance = float(params['cv']['vessel_distance'])
#vessel_diameter = float(params['cv']['vessel_diameter'])
radiusList = [4.6] #mm

for radius in radiusList:
    SphereMesh(radius) 

