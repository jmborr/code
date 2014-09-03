rootdir='/projects/development/VisionElasticFit/tylenol'
bankID='5'
for wid in range(1024):
  swid=str(wid)
  ExtractSingleSpectrum(InputWorkspace='bank'+bankID,WorkspaceIndex=wid,OutputWorkspace='s'+swid)
  Rebin(InputWorkspace='s'+swid,Params='0,4,33600',OutputWorkspace='h'+swid)
  SaveAscii(InputWorkspace='h'+swid,Filename=rootdir+'/histograms/bank'+bankID+'/h'+swid,Separator='Space',CommentIndicator='#')
  DeleteWorkspace('s'+swid)
  DeleteWorkspace('h'+swid)