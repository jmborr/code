'''
This script will try to find the elastic line and the adjacent peak
for a certain number of pixels.
Edit the below variables, then run the script
'''

nexus_file='/projects/development/VisionElasticFit/test2/VIS_1076_event.nxs'  # file containing the run It must be an EVENT nexus file
tofRange=(2000,5000) #range where we expect the elastic line to show up
dt=5.0  # bin spacing (in microseconds)
Max=5000  # maximum time of flight where we expect to find the elastic line
outParmFile='/tmp/collectNevents.dat' # file to output summary fitting data
nbank=20  # number of inelastic banks

############################

from copy import deepcopy

pFile = open(outParmFile,'w')
pFile.write('#1.bank 2.workspaceID 3.Nevents 4.min_nevents\n')
for ibank in (1,2,4,5,7,8,10,11,13,14,16,17,19,20):
  Load(Filename=nexus_file,OutputWorkspace='events',BankName=r'bank%d'%ibank) # load only inelastic banks
  ws=mantid.getMatrixWorkspace('events')
  SumSpectra(InputWorkspace='events',OutputWorkspace='events_summed')
  sws=mantid.getMatrixWorkspace('events_summed')
  min_nevents = sws.readY(0)[0]/ws.getNumberHistograms()  # required minimum number of events  per histogram
  junkLine=''
  for ix in range( ws.getNumberHistograms() ):
    sp = ExtractSingleSpectrum(InputWorkspace='events',WorkspaceIndex=ix,OutputWorkspace='events_single')
    det_id=sp.getDetector(0).getID()
    nevents = sp.readY(0)[0]
    junkLine+='%2d %5d %7d %7d\n'%(ibank,det_id,nevents,min_nevents)
    DeleteWorkspace('events_single')
  DeleteWorkspace('events')
  DeleteWorkspace('events_summed')
  pFile.write(junkLine)
pFile.close()
