from itertools import combinations_with_replacement, product
import numpy as np
import copy


# Given a molecule, find all sets of molecules that could make up the molecule

# Pre: 
# Input allSet, a dictionary containing all molecules to make combination, and molecule, a target to make 
# combination of. 

# Post:
# seperateSmileStrings=False: Returns 1d array of molecular combination jointed with * (suitable for training)
# seperateSmileStrings=True: Returns 2d array of molecular combinations (smile strings are seperated)

def findCombinations(allSet, molecule, seperateSmileStrings=False):
    
    # Helper functions below
    #-------------------------------------------
    
    # Finds the tuple of the given molecule
    # Post: Return the molecule's tuple, return 'Not found' if tuple is not found in given set
    def findTuple():
        for key in allSet:
            arr = allSet.get(key)
            for curr in arr:
                if curr == molecule:
                    return key
        return 'Not found'
    
    def generateTupleCombinations(target, possibleTuples):
        def findCombOf2(targets, options):
            result = []
            combinations = list(combinations_with_replacement(options, 2))
            for comb in combinations:
                first = comb[0]
                second = comb[1]
                total = tuple(map(lambda x, y: x + y, first, second))
                for target in targets:
                    if total == target:
                        tree = {}
                        tree[total] = [first, second]
                        result.append(tree)
            return result

        def generateTargets(trees):
            targets = []
            for tree in trees: # tree is a dic
                key = list(tree.keys())[0]
                leaves = tree[key]
                for leaf in leaves:
                    if leaf not in targets:
                        targets.append(leaf)
            return targets

        # generate next combination set with current combination (Recursive method)
        def generateComb(targetDict, restDict):
            if len(targetDict) != 0:
                result = []
                for comb in targetDict:
                    for molecule in comb:
                        if molecule in restDict:
                            pieces = restDict[molecule]
                            for piece in pieces:
                                newComb = comb.copy()
                                newComb.remove(molecule)
                                newComb += piece
                                newComb.sort()
                                if newComb not in result:
                                    result.append(newComb)
                targetDict += generateComb(result, restDict)
            return targetDict
        
        def collapseTrees(resultTrees):
            result = {}
            for dict in resultTrees:
                curr_key = list(dict.keys())[0]
                if curr_key in result:
                    result[curr_key].append(dict[curr_key])
                else:
                    result[curr_key] = []
                    result[curr_key].append(dict[curr_key])
            return result

        resultTrees = []
        targets = [target]
        while len(targets) != 0:
            trees = findCombOf2(targets, possibleTuples)
            targets = generateTargets(trees)
            for tree in trees:
                resultTrees.append(tree)

        collapedTrees = collapseTrees(resultTrees)
        return generateComb(collapedTrees[target], collapedTrees)
    
    # Convert all tuples into smile strings
    # Pre: Input a 2d list of tuple combinations
    # Post: Output a 2d list of smile combinations
    def tupleCombToSmileComb(combSet):
        result = []
        for tupleComb in combSet:
            smileComb = []
            for i in range(len(tupleComb)):
                smileArr = allSet[tupleComb[i]]
                smileComb.append(smileArr)
            smileComb = product(*smileComb)
            comb = []
            for curr in smileComb:
                result.append(list(curr))
            result.append(comb)
        return result
    
    # Procedures starts below
    # --------------------------------------------
    
    # Search molecule's tuple
    atomCount = findTuple()
    if atomCount == 'Not found':
        raise Exception('Molecule not found in given database')
    possibleTuples = list(allSet.keys())
    
    # Remove tuples with empty value
    remove = []
    for curr in possibleTuples:
        if allSet[curr] == ['']:
            remove.append(curr)
            
    # Put all tuple that have C, O or H exceeding the given tuple in remove list
    for key in possibleTuples:
        for i in range(len(key)):
            if key[i] > atomCount[i]:
                remove.append(key)
                
    # Put isomers remove list, no need to make combinations with them. Add them back in later in result
    remove.append(atomCount)
    
    # Actual removing step
    for curr in remove:
        if curr in possibleTuples:
            possibleTuples.remove(curr)

    
    combSet = generateTupleCombinations(atomCount, possibleTuples)
    
    
    # Convert tuple to smile strings
    smileCombSet = tupleCombToSmileComb(combSet)
    smileCombSet = list(filter(lambda a: a != [], smileCombSet)) # Filter out werid empty list made in product()
    # Add isomers into smileCombSet
    for smile in allSet[atomCount]:
        smileCombSet.append([smile])
    
    if seperateSmileStrings == True:
        return smileCombSet
    else: # Join them with *
        result = []       
        for comb in smileCombSet:
            combString = comb[0]
            for i in range(1,len(comb)):
                combString += '.' + comb[i]
            result.append(combString)
        return result