energy_keys='ETITLE: TS BOND ANGLE DIHED IMPRP ELECT VDW BOUNDARY MISC KINETIC TOTAL TEMP POTENTIAL TOTAL3 TEMPAVG PRESSURE GPRESSURE VOLUME PRESSAVG GPRESSAVG'.split()

def fetch_energy_column(logfile,column):
    """return the colum in the ETITLE: lines

    Scan the log file for each line beginning with 'ENERGY:', split it, and save the requested column index
    
    Args:
      logfile (string): log file pathname
      column (int): colum index in the ETITLE: line, splitted with split()

    Returns: 
      the column in ETITLE: as a list

    Raises:
      n/a
    """
    print 'Fetch values for ',energy_keys[column]
    values=[]
    for line in open(logfile).readlines():
        if line[:7]=='ENERGY:': 
            fields=line.split()
            if column < len(fields):
                values.append(float(fields[column]))
            else:
                values.append(0)
    return values

def plot_energy_column(logfile,column):
    """Simple plot of a column in ETITLE: lines

    Args:
      logfile (string): log file pathname
      column (int): colum index in the ETITLE: line, splitted with split()

    Returns: 
      plot object

    Raises:
      n/a
    """
    import matplotlib
    matplotlib.pyplot.xlabel('time')
    matplotlib.pyplot.ylabel(energy_keys[column])
    x = matplotlib.pyplot.plot( fetch_energy_column(logfile,1), fetch_energy_column(logfile,column))
    matplotlib.pyplot.show()
    return x
    
