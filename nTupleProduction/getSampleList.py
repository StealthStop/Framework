import dict_master
from optparse import OptionParser

def main():

	parser = OptionParser()
	parser.add_option("-y", "--year", dest="year", default = "Autumn18", action="store", help="Indicate the desired year. Options: Summer16v3, Fall17, Autumn18")
	parser.add_option("-s", "--sampleType", dest="sampType", default = "QCD_HT", action="store", help="Indicate the sample type desired. Options: QCD_HT, QCD_PT, QCD_PT_MuEnriched, diboson, dyjets, gjets, ttbar, tth, wjets, zjets, singlet")
	parser.add_option("-a", "--allSamps", dest = "allSamps", default = False, action="store_true", help="List all sample types for specified year")	
	parser.add_option("-p", "--print", dest = "doPrint", default = False, action="store_true", help="Print sample names in terminal")
        parser.add_option("-f", "--saveFile", dest = "saveFile", default = False, action="store_true", help="Save sample names to .txt file")
	(options, args) = parser.parse_args()

	year = options.year
	sampType = options.sampType
	allSamps = options.allSamps
	doPrint = options.doPrint
        saveFile = options.saveFile
	listToPrint = []
	if allSamps:
		for s in dict_master.flist[year]:
			listToPrint.extend([i[0] for i in dict_master.flist[year][s]])
		if doPrint:
			for l in listToPrint:
				print("%s\n" % l)
                if saveFile:
                        with open("sampleList_"+year+"_All.txt", "w") as f:
                                for l in listToPrint:
                                        f.write("%s\n" % l)
                        f.close()
                        print("\n" "Saving all sample names for "+year+" to sampleList_"+year+"_All.txt" "\n")
	else:
		listToPrint = [i[0] for i in dict_master.flist[year][sampType]]
		if doPrint:
			for l in listToPrint:
		 		print("%s\n" % l)
		if saveFile:
                        with open("sampleList_"+year+"_"+sampType+".txt", "w") as f:
                                for l in listToPrint:
                                        f.write("%s\n" % l)
			f.close()
                        print("Saving all sample names for "+year+" to sampleList_"+year+"_"+sampType+".txt" "\n")
	

if __name__ == "__main__":
	main()
