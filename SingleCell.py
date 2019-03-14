# Import Packages
import pandas as pd
import numpy as np



class SingleCell:
    
    dataset = None
    data = None
    celldata = None
    genedata = None
    dim = None

    excl_cells = None # Boolean vector containing cells that are excluded from analysis (0-excluded cells)
    excl_genes = None # Boolean vector containing genes that are excluded from analysis (0-excluded genes)
    spike_ins_names = np.array([])
    
    def __init__(self, dataset, data, celldata, genedata):
        
        self.dataset = dataset
        self.data = data
        self.celldata = celldata
        self.genedata = genedata
        self.dim = data.shape
        
        self.excl_cells = np.ones(self.dim[1], dtype=bool)
        self.excl_genes = np.ones(self.dim[0], dtype=bool)

    def getLogCounts(self):
        data = self.data.values
        data = data[self.excl_genes, :]
        data = data[:, self.excl_cells]
        return np.log2(data + 1)
    
    def getCounts(self):
        data = self.data.values
        data = data[self.excl_genes, :]
        data = data[:, self.excl_cells]
        return data
    
    def setExcludedCells(self, cells, column):
        # Exclude the following cells
        cell_labels = self.celldata[column].values
        for i in cells:
            for j in range(cell_labels.size):
                if (i == cell_labels[j]):
                    self.excl_cells[j] = False
                    
                    
    def setExcludedGenes(self, genes):
        # Exclude the following genes
        gene_labels = self.genedata["feature_symbol"].values
        for i in genes:
            for j in range(gene_labels.size):
                if (i == gene_labels[j]):
                    self.excl_genes[j] = False
                    
            
    
    def checkCellDataColumn(self, column):
        cell_data_cols = self.celldata.columns.values
        
        status = False
        
        for i in cell_data_cols:
            if (i == column):
                status = True
                
        return status
    

    def checkGeneDataColumn(self, column):
        gene_data_cols = self.genedata.columns.values

        status = False

        for i in gene_data_cols:
            if (i == column):
                status = True

        return status


    def isSpike(self, spike_type):
        
        mask = np.zeros(self.dim[0], dtype=bool)
        
        if (self.checkGeneDataColumn("feature_symbol")):
            gene_labels = self.genedata["feature_symbol"].values
            mask = pd.Series(gene_labels).str.contains(spike_type).tolist()
            
            # ser = self.genedata.duplicated(subset="feature_symbol")
            total = np.sum(mask)
            if (total > 0):
                print("Spike-ins '", spike_type, "' exists in the dataset!")
                self.spike_ins_names = np.append(self.spike_ins_names, spike_type)
                self.excl_genes[mask] = False
            else:
                print("Spike-ins '", spike_type, "' does not exist in the dataset!")
            
        else:
            print("'feature_symbol' column not found in gene data!")
        
    
        
    def addGeneDataColumn(self, col_data, col_name):
        self.removeGeneDataColumn(col_name)
        
        df = pd.DataFrame(data = col_data, columns = [col_name], index = self.genedata.index)
        self.genedata = pd.concat([self.genedata, df], axis = 1)
        


    def addCellDataColumn(self, col_data, col_name):
        self.removeCellDataColumn(col_name)
        
        data = np.zeros(self.dim[1])
        data[self.excl_cells] = col_data
        
        df = pd.DataFrame(data = data, columns = [col_name], index = self.celldata.index)
        self.celldata = pd.concat([self.celldata, df], axis = 1)
        


    def removeCellDataColumn(self, column):
        
        if (self.checkCellDataColumn(column)):
            self.celldata = self.celldata.drop([column], axis=1)
            print ("Removing '", column, "' from CellData assay")

        else:
            print ("'", column, "' does not exist in CellData assay")

    

    def removeGeneDataColumn(self, column):
        
        if (self.checkGeneDataColumn(column)):
            self.genedata = self.genedata.drop([column], axis=1)
            print ("Removing '", column, "' from GeneData assay")
    
        else:
            print ("'", column, "' does not exist in GeneData assay")


    """
    ToDo: Check if column exists in the relevant assay
    """
    def getCellColumnValues(self, column):
        
        return self.celldata[column].values[self.excl_cells]
    
    def getGeneColumnValues(self, column):
        
        return self.genedata[column].values[self.excl_genes]
           

  
    def getDistinctCellTypes(self, column):
        
        cellarray = self.getCellColumnValues(column)
        cells = cellarray

        cell_types = np.array([]) # Distinct cell labels in the cellarray

        # Find distinct cell lables in the cell array
        while (cells.size != 0):
            cell_types = np.append(cell_types, cells[0])
            mask = cells != cells[0]
            cells = cells[mask]

        return cell_types



    def getNumericCellLabels(self, column):
        
        cellarray = self.getCellColumnValues(column)
        cell_types = self.getDistinctCellTypes(column)

        num_cell_labels = np.zeros(cellarray.size) # Numeric cell labels

        # Find numeric cell labels for the cell array 
        for i in range(cell_types.size):
            mask = cellarray == cell_types[i]
            num_cell_labels[mask] = i
        
        return num_cell_labels
        
        


