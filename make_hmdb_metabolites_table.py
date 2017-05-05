# Adam Berman
# Princeton University
# 5 May 2017


import os
import xmltodict
import xml.etree.ElementTree as ET
import itertools


# Parse metabolites

all_metabolites = []
all_inchi_metabolites = []
all_blank_metabolites = []


# Iterate over all metabolites
# Folder containing all metabolite's xml files (named "hmdb_metabolites") can be found 
# at the following link: http://www.hmdb.ca/system/downloads/current/hmdb_metabolites.zip
for filename in os.listdir('/Users/adamberman/Independent_Work/hmdb_metabolites'):

	if filename.startswith('HMDB'):

		# Parse the metabolite
		root = ET.parse('/Users/adamberman/Independent_Work/hmdb_metabolites/' + filename).getroot()
		origin = root.find('ontology').find('origins').find('origin')

		# Keep only endogenous metabolites
		if (origin != None) and (origin.text == 'Endogenous'):

			# Collect all genes associated with the given metabolite
			gene_names = []
			for protein in root.find('protein_associations').findall('protein'):
				gene_name = protein.find('gene_name').text
				if gene_name != None:
					gene_names.append(gene_name.strip())

			# Only continue if the metabolite has at least one associated
			if gene_names:

				# Collect the name and KEGG ID of the metabolite
				name = root.find('name').text.strip()
				kegg = root.find('kegg_id').text
				if kegg == None:
					kegg = "None"
				else:
					kegg = kegg.strip()

				# Collect the synonyms of the metabolite
				synonyms = []
				for synonym in root.find('synonyms').findall('synonym'):
					synonyms.append(synonym.text.strip())

				# Collect the SMILES identifier of the metabolite
				smiles = root.find('smiles').text
				if smiles != None:
					smiles = smiles.strip()
				
				# InChI instead of SMILES
				if smiles == None:
					
					inchi = root.find('inchi').text
					if inchi != None:
						inchi = inchi.strip()

					if inchi != None: 
						inchi_key = root.find('inchikey').text
						if inchi_key != None:
							inchi_key = inchi_key.strip()

						# Add new metabolite to all_inchi_metabolites
						metabolite = {'name': name, 'synonyms': synonyms, 'kegg': kegg, 'gene_names': gene_names, 'inchi': inchi, 'inchi_key': inchi_key, 'filename': filename}
						all_inchi_metabolites.append(metabolite)

					else:
						# Add new metabolite to all_blank_metabolites
						metabolite = {'name': name, 'synonyms': synonyms, 'kegg': kegg, 'gene_names': gene_names, 'filename': filename}
						all_blank_metabolites.append(metabolite)


				# SMILES
				else:
					# Add new metabolite to all_metabolites
					metabolite = {'name': name, 'synonyms': synonyms, 'kegg': kegg, 'gene_names': gene_names, 'smiles': smiles, 'filename': filename}
					all_metabolites.append(metabolite)



# Remove synonymous metabolites from all_metabolites
smiles_dictionary = {}
for metabolite in all_metabolites:

	genes = metabolite.get('gene_names')
	smls = metabolite.get('smiles')
	kgg = metabolite.get('kegg')
	syn = metabolite.get('synonyms')
	name = metabolite.get('name')
	name_and_syn = syn + [name]

	merged = False
	for key, value in smiles_dictionary.iteritems():

		# Merge metabolite with existing smiles_dictionary entry if 
		# synonyms overlap, or KEDD ID's are the same, or SMILES are the same
		# or ((len(genes) > 1) and (genes == value.get('gene_names')))
		if ((key in name_and_syn) 
			or (not set(value.get('synonyms')).isdisjoint(name_and_syn)) 
			or ((kgg != "None") and (kgg == value.get('kegg'))) 
			or (smls == value.get('smiles'))):
			
			# Perform merge operation
			new_syn = list(set(name_and_syn + value.get('synonyms')))
			if key in new_syn:
				new_syn.remove(key)

			new_gn = list(set(metabolite.get('gene_names') + value.get('gene_names')))
			new_gn.sort()

			# Create new entry and store it into the dictionary
			new_value = {'synonyms': new_syn, 'kegg': value.get('kegg'), 'gene_names': new_gn, 
							'smiles': value.get('smiles'), 'filename': value.get('filename')}

			smiles_dictionary[key] = new_value

			# Denote that merge occured
			merged = True
			break

	# If the metabolite was not merged, create new smiles_dictionary entry for metabolite
	if not merged:
		smiles_dictionary[name] = {'synonyms': syn, 'kegg': metabolite.get('kegg'), 
									'gene_names': metabolite.get('gene_names'), 
									'smiles': metabolite.get('smiles'), 
									'filename': metabolite.get('filename')}



# Print SMILES metabolites into table
print 'Metabolite' + '\t' + 'SMILES' + '\t' + 'KEGG' + '\t' + 'Gene Name'
for key, value in smiles_dictionary.iteritems():
	n = key
	s = value.get('smiles')
	k = value.get('kegg')
	gns = value.get('gene_names')

	for gn in gns:
		#print '{0:40s} {1:500s} {2:30s}'.format(n, s, gn)
		print n + '\t' + s + '\t' + k + '\t' + gn






# TEST CODE

'''
# Simple test version of all_metabolites
all_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'hound', 'synonyms': ['new', 'newer', 'newest'], 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'gene_names': ['GGG', 'HHH'], 'smiles': 'IOU1K9_4', 'filename': 'Dogger_4.xml'}]
'''

'''
# Simple test version of all_inchi_metabolites
all_inchi_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'gene_names': ['1AA', 'B1B'], 'inchi': 'inchi_1', 'inchi_key': 'inchi_key_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'gene_names': ['B2B', 'BBB'], 'inchi': 'inchi_2', 'inchi_key': 'inchi_key_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'hound', 'synonyms': ['new', 'newer', 'newest'], 'gene_names': ['EEE', 'FFF'], 'inchi': 'inchi_3', 'inchi_key': 'inchi_key_3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'gene_names': ['GGG', 'HHH'], 'inchi': 'inchi_4', 'inchi_key': 'inchi_key_4', 'filename': 'Dogger_4.xml'}]
'''

'''
all_metabolites = [{'name': 'dog', 'synonyms': ['canine', 'pupper', 'poochy'], 'kegg': 'dog_kegg', 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'hound', 'synonyms': ['dog', 'pupper', 'doggo'], 'kegg': 'hound_kegg', 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'Title', 'synonyms': ['new', 'newer', 'newest'], 'kegg': 'Title_kegg', 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'filename': 'Dogger_3.xml'},
					{'name': 'New', 'synonyms': ['Title', 'hound', 'poochy'], 'kegg': 'New_kegg', 'gene_names': ['GGG', 'HHH'], 'smiles': 'IOU1K9_4', 'filename': 'Dogger_4.xml'}]
'''


'''
all_metabolites = [{'name': 'test1', 'synonyms': ['a', 'b', 'c'], 'kegg': 'same_kegg', 'gene_names': ['1AA', 'B1B'], 'smiles': 'IOU1K9_1', 'filename': 'Dogger_1.xml'}, 
					{'name': 'test2', 'synonyms': ['d', 'e', 'f'], 'kegg': 'diff_kegg', 'gene_names': ['B2B', 'BBB'], 'smiles': 'IOU1K9_2', 'filename': 'Dogger_2.xml'}, 
					{'name': 'test3', 'synonyms': ['g', 'h', 'i'], 'kegg': 'same_kegg', 'gene_names': ['EEE', 'FFF'], 'smiles': 'IOU1K9_3', 'filename': 'Dogger_3.xml'},
					{'name': 'test4', 'synonyms': ['j', 'k', 'l'], 'kegg': 'same_kegg', 'gene_names': ['GGG', 'HHH'], 'smiles': 'IOU1K9_4', 'filename': 'Dogger_4.xml'}]
'''
