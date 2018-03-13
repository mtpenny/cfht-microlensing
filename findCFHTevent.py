import re
import numpy as np

eventfiles = {'O':'ews_cfht.txt', 'M':'moa_cfht.txt', 'K':'kmt_cfht.txt'}

#Function for loading event data from file
def event_from_file(event,filename):
    
    elf = open(filename,'r')
    regex = re.compile('^%s ' % (event))
    ret = []
    for line in elf:
        if regex.match(line):
            line = re.sub('Cccd',' ccd',line) #Split the field and chip code
            ret = line.split()
            break
    elf.close()
    return ret


def find_event(event,field=None):

    '''Search for an event in an event list file (format event ra dec 
       field-chip -code) and parse the code. If field specified, requires
       that the event be in that field.

       Return: A dictionary containing the matched event name (event), 
               ra, dec, field and chip (if available).'''

    dictret={'event':None, 'ra':None, 'dec':None, 'field':None, 'chip':None}

    #Make sure MOA events have the right format
    if event[0] == 'M' and event[4] != 'M':
        event = event[0:4] + 'M' + event[4:]
    
    ret=[]
    if event[0] in eventfiles:
        #We know where the event should be
        ret = event_from_file(event,eventfiles[event[0]])
    else:
        print event, field
        for key,val in eventfiles:
            tmpret = event_from_file(event,val)
            if tmpeventdata!=None:
                ret=tmpret
                break

    #Find the event info in the event file

    #Deal with the case where event wasn't found (user handles it)
    if len(ret)==0:
        #Event wasn't found, return dict of Nones
        return dictret
    
    if len(ret)<=3:
        #Event was found, but isn't on a chip
        dictret['event']=ret[0]
        dictret['ra']=ret[1]
        dictret['dec']=ret[2]
        return dictret

    #Deal with the case where event has multiple fields
    if len(ret)>5:
        #Event landed in more than one field
        if field != None:
            #User asked for a field, search for it, otherwise we'll just use
            #the first
            for i,tmpfield in enumerate(ret[3::2]):
                if tmpfield == field:
                    ret[3]=tmpfield
                    ret[4]=ret[4+2*i]

    if field != None and ret[3] != field:
        #Event is not in specified field
        dictret['event']=ret[0]
        return dictret

    dictret['event']=ret[0]
    dictret['ra']=np.float(ret[1])
    dictret['dec']=np.float(ret[2])
    dictret['field']=ret[3]
    dictret['chip']=ret[4]
    
    return dictret

def find_event_errors(eventdata):
    if not 'event' in eventdata:
        return "eventdata not initialized by call to find_event"
    if eventdata['event']==None:
        return "Event not found"
    if eventdata['field']==None:
        if eventdata['ra']==None:
            return "Event %s is not in the specified field" % (eventdata['event'])
        else:
            return "Event %s does not land on a chip" % (eventdata['event'])

    return None
    
