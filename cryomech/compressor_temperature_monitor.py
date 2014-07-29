#
#   Run this program
#
#   Don't try to mess with it
#
#   It's really complicated
#
#
#     Rob H
#

import binascii
import time
from visa import *
import csv

output_filename = 'comp_temps-'+str(int(time.time()))+'.txt'

#---------------------

Comp_port = 'ASRL3'
Comp = instrument(Comp_port, baud_rate=115200, term_chars = CR, values_format = single | big_endian)

addr = '\x10'
stx = '\x02'
cmd_rsp = '\x80'
cr = '\x0D'

waittime = 2

read = '\x63'
write = '\x61'

hashcode = {
'CODE_SUM': '\x2B\x0D',
'MEM_LOSS': '\x80\x1A',
'CPU_TEMP': '\x35\x74',
'BATT_OK': '\xA3\x7A',
'BATT_LOW': '\x0B\x8B',
'COMP_MINUTES': '\x45\x4C',
'MOTOR_CURR_A': '\x63\x8B',
'RI_RMT_COMP_START': '\xBA\xF7',
'RI_RMT_COMP_STOP': '\x3D\x85',
'RI_RMT_COMP_ILOK': '\xB1\x5A',
'RI_SLVL': '\x95\xE3',
'TEMP_TNTH_DEG': '\x0D\x8F',
'TEMP_TNTH_DEG_MINS': '\x6E\x58',
'TEMP_TNTH_DEG_MAXES': '\x8A\x1C',
'CLR_TEMP_PRES_MMMARKERS': '\xD3\xDB',
'TEMP_ERR_ANY': '\x6E\x2D',
'PRES_TNTH_PSI': '\xAA\x50',
'PRES_TNTH_PSI_MINS': '\x5E\x0B',
'PRES_TNTH_PSI_MAXES': '\x7A\x62',
'CLR_TEMP_PRES_MMMARKERS': '\xD3\xDB',
'PRES_ERR_ANY': '\xF8\x2B',
'H_ALP': '\xBB\x94',
'H_AHP': '\x7E\x90',
'H_ADP': '\x31\x9C',
'H_DPAC': '\x66\xFA',
'DIODES_UV': '\x8E\xEA',
'DIODES_TEMP_CDK': '\x58\x13',
'DIODES_ERR': '\xD6\x44',
'DCAL_SEL': '\x99\x65',
'EV_START_COMP_REM': '\xD5\x01',
'EV_STOP_COMP_REM': '\xC5\x98',
'COMP_ON': '\x5F\x95',
'ERR_CODE_STATUS': '\x65\xA4'
}

escapes = {
	'\x02': '\x07\x30',
	'\x0D': '\x07\x31',
	'\x07': '\x07\x32'
}
#first thing: datastructure
#
#STX	ADDR	CMD_RSP		[DATA...]	CKSUM1	CKSUM2	CR
#
#STX 		\x02
#ADDR		address of the slave to be spoken to?
#CMD_RSP	send \x80
#DATA 		must be escaped: 
#			\x02 -> \x07\x30 
#			\x0D -> \x07\x31
#			\x07 -> \x07\x32
#CKSUM1		compute mod-256 checksum of ADDR, CMD_RSP, and DATA
#			CKSUM1 is upper 4 bits 'read as a nibble' + \x30, w/ result \x30 to \x3F
#CKSUM2		is lower 4 bits + \x30 w/ result \x30 to \3F
#CR: 		\x13 or '\r' (I think)

def escape_data(data,escapes):
	result = ''
	for c in data:
		if c in escapes:
			result+=escapes[c]
		else:
			result+=c
	return result

def capture_data(result,escapes):
        '''
        does the reverse of escape_data - recombines it
        '''
        stx = result[0]
        addr = result[1]
        cmd_rsp = result[2]
        scksum1 = result[-2]
        scksum2 = result[-1]
        data = result[3:-2]
        for key in escapes.keys():
                data = data.replace(escapes[key],key)
        chk1,chk2 = checksum(addr,cmd_rsp,data)
        if scksum1!=chk1 or scksum2 != chk2:
                print 'checksum failure, bad data acquisition'
                return False
        else:
                return int(binascii.hexlify(data[4:]),16)*0.1

def gen_data(readwrite,hashcode,index='\x00',writedata=''):
	if readwrite == '\x63': outlen = 4
	if readwrite == '\x61': outlen = 8
	output = readwrite + hashcode + index + writedata
	if len(output) == outlen:
		return output
	else:
		print 'readwrite suggests length',outlen,'but gen_data input only',len(output),':',binascii.hexlify(output)


def checksum(addr, cmd_rsp, data):
	result = 0
	for c in addr + cmd_rsp + data:
		result += int(binascii.hexlify(c),16)
		#print binascii.hexlify(c),hex(result)
	result = binascii.hexlify(chr(result%256))
	#print result
	chk1 = binascii.unhexlify('3'+ result[0])
	chk2 = binascii.unhexlify('3'+ result[1])
	return chk1, chk2


def sendstring(stx,addr,cmd_rsp,escdata,cksum1,cksum2,cr):
	result = stx+addr+cmd_rsp+escdata+cksum1+cksum2+cr
	return result


def gen_command(name,hashcode,escapes,index='\x00',writedata='',readwrite='\x63'):
        data = gen_data(readwrite,hashcode[name],index=index,writedata=writedata)
        escdata = escape_data(data,escapes)
        chk1, chk2 = checksum(addr,cmd_rsp,data)
        command = sendstring(stx,addr,cmd_rsp,escdata,chk1,chk2,cr)
        return command


name = 'TEMP_TNTH_DEG'
indexes = ['\x00','\x01','\x02','\x03']
#writedata = '\x00\x00\x00\x01' #only used when readwrite = write
writedata = ''
readwrite = read

#
# generate your command using:
#
# gen_command(name_of_function, hashcode_dictionary, escapecodes_dictionary,
#       and optionally the index, data to be written if writing not reading, and
#       whether you're performing a read or a write)
#
# this generates the packet to be sent to the compressor via RS232
#
# once you have this you can use the PyVISA command comp.ask() to interface
# with the compressor
#

fileHandle = open(output_filename,'a')
csvHandle = csv.writer(fileHandle,delimiter=',')
csvHandle.writerow(['time (s)','asctime (human readable)',
                    'input water temp (C)','output water temp (C)','helium temp (C)','oil temp (C)'])
fileHandle.close()

commands = []
for index in indexes:
        commands.append(gen_command(name,hashcode,escapes,index=index,writedata=writedata,readwrite=readwrite))

try:
        while True:
                results = []
                results.append(time.time())
                results.append(time.asctime())
                for command in commands:
                        while True:
                                response = Comp.ask(command)
                                captured = capture_data(response,escapes)
                                if captured != False:
                                        results.append(captured)
                                        break
                for cell in results[1:]:
                        print cell,'\t',
                print '\n',
                fileHandle = open(output_filename,'a')
                csvHandle = csv.writer(fileHandle,delimiter=',')
                csvHandle.writerow(results)
                fileHandle.close()
                time.sleep(waittime)
except KeyboardInterrupt:
        print '\nClosing gracefully...'
        try:
                fileHandle.close()
        except:
                pass


