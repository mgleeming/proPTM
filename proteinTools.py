import os, sys, re
from pyteomics import mass, parser

class Modification (object):
    def __init__(self, name,  netFormulaChange = None, reactiveResidues = None):
        self.name = name
        self.netFormulaChange = netFormulaChange
        self.netMassChange = mass.calculate_mass(formula = netFormulaChange)
        self.reactiveResidues = reactiveResidues
        return
    def __repr__(self):
        return (4*'%s, ' %(self.name, self.netFormulaChange, self.netMassChange, self.reactiveResidues))


class Protein (object):
    def __init__(self, sequence = None, PDB = None, minLength = 6):
        if sequence:
            self.sequence = ''.join(re.split("[^a-zA-Z]*", sequence))
        if PDB:
            pass # get PDB file
        self.minLength = minLength
        self.peptides = []
        self.modifications = []
        return

    def digest(self, missedCleavages = 0, charges = [1,2,3]):
        if self.sequence:
            fragments = parser.cleave(self.sequence, parser.expasy_rules['trypsin'], missedCleavages)
            for fragment in fragments:
                if len(fragment) < self.minLength:
                    continue
                self.peptides.append(
                    Peptide (fragment, self.modifications, charges)
                )
        return

    def applyModification(self, newMod):
        if isinstance(newMod, Modification):
            self.modifications.append(newMod)
        else:
            print 'First create Modification instance'
        return

    def getModifications(self):
        for i in range(len(self.modifications)):
            print self.modifications[i]
        return



class Peptide (object):
    def __init__(self, sequence, proteinModifications, charges):
        self.sequence = sequence
        self.proteinModifications = proteinModifications
        self.charges = charges
        self.nativeMass = self.getNativeMass()
        self.modifiedMass = self.getModifiedMass()
        #self.fragments = self.getPeptideFragments()
        self.pseudoMolecularIons = self.getPseudoMolecularIons()
        return

    def getNativeMass(self):
        return mass.fast_mass(sequence=self.sequence)

    def getModifiedMass(self):
        # get list of modifications in peptides
        totalModsMass = 0
        for i in range(len(self.sequence)):
            for m in self.proteinModifications:
                if self.sequence[i].upper() in m.reactiveResidues:
                    totalModsMass += m.netMassChange
        return totalModsMass + self.nativeMass

    def getPseudoMolecularIons(self):
        ions = {}
        for i in self.charges:
            ions[i] = (self.modifiedMass + i * 1.00728) / float(i)
        return ions



'''

DCV = pptm.Modification(netFormulaChange = 'C2O3PH5', reactiveResidues = ['S', 'T', 'Y'])

HSA = pptm.Protein(sequence = 'PEPTIDE')
HSa.applyModification(DCV)



'''

