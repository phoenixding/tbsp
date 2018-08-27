import pdb,sys,os

class BioList(list):
	def ex2File(self,outname):
                if len(self)>0 and (type(self[0])==list):
                    self=[[str(item) for item in k] for k in self]
                    self=["\t".join(item) for item in self]
                out=[str(item) for item in self]
                out='\n'.join(out)
                f=open(outname,'w')
                f.write(out)
                f.close()
		
		
class bed:
	def __init__(self,bedRows):		
		self.bedType=len(bedRows[0])
		self.bedRows=bedRows
	
	def bed2Wig(self):
		# use with caution: all bed lengths must be for same. 
		if self.bedType<5:
			raise ValueError("requires >bed5")
		#
		print("building dictionary")
		dchr={}
		for row in self.bedRows:
			chr_row=row[0]
			start_row=row[1]
			end_row=row[2]
			name_row=row[3]
			score_row=row[4]
			if chr_row not in dchr:
				dchr[chr_row]=[[start_row,end_row,score_row]]
			else:
				dchr[chr_row].append([start_row,end_row,score_row])
		
		#
		span=int(self.bedRows[0][2])-int(self.bedRows[0][1])
		print("exporting wig")
		outstr=""
		for chri in dchr:
			outstr+="variableStep chrom=%s span=%s"%(chri,span)+"\n"
			chrValues=dchr[chri]
			#pdb.set_trace()
			for j in chrValues:
				#pdb.set_trace()
				jstart=j[0]
				jend=j[1]
				jscore=j[2]
				outstr+=str(jstart)+"\t"+str(jscore)+"\n"
			print(chri)
		return outstr
				
			
