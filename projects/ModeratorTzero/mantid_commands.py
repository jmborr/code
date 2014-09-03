LoadEventNexus(Filename='/SNS/VIS/2012_1_16B_SCI/1/1085/NeXus/VIS_1085_event.nxs',OutputWorkspace='ws',BankName='bank4')

#Apply ModeratorTzero to the events
ModeratorTzero(InputWorkspace='ws',OutputWorkspace='wsx')

#Apply ModeratorTzero to the histograms
Rebin(InputWorkspace='ws',Params=[0,4,16000],OutputWorkspace='hws',PreserveEvents=0)
ModeratorTzero(InputWorkspace='hws',OutputWorkspace='hwsx')

#compare original with transformed events  and with transformed histograms
SumSpectra(InputWorkspace='ws',OutputWorkspace='sum')
SumSpectra(InputWorkspace='wsx',OutputWorkspace='sumx')
SumSpectra(InputWorkspace='hwsx',OutputWorkspace='transformed_histograms')
Rebin(InputWorkspace='sum',Params=[0,4,16000],OutputWorkspace='original')
Rebin(InputWorkspace='sumx',Params=[0,4,16000],OutputWorkspace='transformed_events')
plot((mtd['original'],mtd['transformed_events'],mtd['transformed_histograms']),0)