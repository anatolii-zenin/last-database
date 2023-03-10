import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import re
from glob import glob
import mysql.connector

class DatabaseHandler:
    def __init__(self, folder = "", tablename = ""):
        self.folder = folder
        if self.folder != "":
            self.df = self.readFolder(folder)
            self.writeDatabase(self.df, tablename)

    # setters/getters
    def getFolder(self):
        return self.folder

    def setFolder(self, name):
        self.folder = name

    def getColNames(self, table):
        if table == "poldata":
            return ['date', 'time', 'MJD', 'source', 'RA', 'DEC', 'position',
                           'x', 'y', 'set', 'star', 'origin', 'Q',
                           'Qerr', 'U', 'Uerr', 'PD', 'PDerr', 'PA', 'PAerr']
        elif table == "diffphot":
            return ['date', 'time', 'exptime', 'name', 'aperture_set',
                                  'letters', 'sky', 'fwhm', 'ellipticity', 'value2',
                                  'value3', 'mag', 'mag_error', 'status', 'd_object',
                                  'd_grid_pos', 'd_set', 'd_pol']

    # reads all .instr files in a folder
    # and combines them in a single dataframe
    def readFolder(self, folder = ""):
        if folder == "":
            folder = self.folder
        self.setFolder(folder)
        files = glob(folder + "*.instr")
        print(files)
        dfnames = []
        for i in range(len(files)):
            dfnames.append(str(i))
        dfdict = {name: pd.DataFrame() for name in dfnames}
        dflist = []
        for i in range(len(files)):
            dfdict[dfnames[i]] = self.loadData(files[i])
            dfdict[dfnames[i]]['d_set'] = dfdict[dfnames[i]]['d_set'].astype(int)
            dfdict[dfnames[i]]['name'] = dfdict[dfnames[i]]['name'].astype(str)
            dflist.append(dfdict[dfnames[i]])
            print(files[i], " read, ", len(dfdict[dfnames[i]]), " added" )
        df = pd.concat(dflist, axis=0)
        df.reset_index(drop=True, inplace=True)
        return df

    # auxiliary method for reading files
    def loadData(self, input_file):
        file_name=input_file
        data1=pd.read_csv(input_file, skiprows=lambda x: x%2 == 1 or x < 2, header=None, sep='\s+')
        data2=pd.read_csv(input_file, skiprows=lambda x: x%2 == 0 or x < 2, header=None, sep='\s+')
        data=pd.concat([data1,data2], axis=1)
        data.columns=['date','time','exptime','name','aperture_set','letters',
                       'sky','fwhm','ellipticity','value2','value3','mag',
                       'mag_error','status']
        n_names = len(list(data['name'].str.split('_', 1, expand=True)))
        d_object = data['name'].str.rsplit('_', 1, expand=True)
        params = d_object[1]
        data['d_object']='nan'
        data['d_grid_pos']='nan'
        data['d_set']='nan'
        data['d_pol']='nan'
        for i in range(0,len(d_object)):
            par_i = params[i].split('-')
            data.loc[i, 'd_object'] = d_object[0][i]
            if n_names==2:
                data.loc[i, 'd_grid_pos'] = par_i[0]
                data.loc[i, 'd_set'] = par_i[1]
                data.loc[i, 'd_pol'] = par_i[2]
            else:
                exit()
        data['d_pol']=data['d_pol'].str.replace(".fit",'', regex=True)
        pols = np.array(data['d_pol'], dtype=str)
        pols_num = np.zeros(len(pols), dtype = int)
        for i in range(len(pols)):
            pols_num[i] = int(re.findall(r'\d+', pols[i])[0])
        data['d_pol']=pols_num
        return data

    # methods for calculations
    def getQU(self, I0, I90, I0err, I90err):
        Q = (I0 - I90)/(I0 + I90)
        Qerr = 1/(I0+I90)**2 * np.sqrt((2*I0)**2 * I90err**2 +
                            (2*I90)**2 * I0err**2)
        return (Q, Qerr)

    def getPD(self, Q, U, Qerr, Uerr):
        PD = np.sqrt(Q**2 + U**2)
        PDerr = np.sqrt((Q/np.sqrt(Q**2 + U**2))**2 * Qerr**2 +
                                  (U/np.sqrt(Q**2 + U**2))**2 * Uerr**2)
        return (PD, PDerr)

    def getPA(self, Q, U, Qerr, Uerr):
        PA = np.rad2deg(np.arctan(U/Q))
        PAerr = np.rad2deg(np.sqrt((U/(Q**2 + U**2))**2 * Qerr**2 +
                                  (Q/(Q**2 + U**2))**2 * Uerr**2))
        return (PA, PAerr)

    def getPol(self, I0, I45, I90, I135, I0err, I45err, I90err, I135err):
        Q, Qerr = self.getQU(I0, I90, I0err, I90err)
        U, Uerr = self.getQU(I45, I135, I45err, I135err)
        PD, PDerr = self.getPD(Q, U, Qerr, Uerr)
        PA, PAerr = self.getPA(Q, U, Qerr, Uerr)
        return PD, PDerr, PA, PAerr, Q, Qerr, U, Uerr

    def updatePolTable(self, tablename_in = "diffphot", tablename_out = "poldata"):
        if tablename_in != "diffphot":
            raise Exception("makePolDataframe: unknown input type: %s" % tablename_in)
            exit()
        if tablename_out != "poldata":
            raise Exception("makePolDataframe: unknown output type: %s" % tablename_out)
            exit()
        df = self.readDatabase(tablename_in)
        df_pol = self.calcPol(df)
        self.writeDatabase(df_pol, tablename_out)


    # make a polarisation dataframe using input dataframe
    def calcPol(self, source, aperture = 'OBJECT_O_BAD'):
        cols = self.getColNames("poldata")
        df_pol = pd.DataFrame(columns = cols)
        I = {}
        I_err = {}
        dates = source['date'].unique()
        # loops through all the parameters:
        # 1. date
        # 2. objects
        # 3. star
        # 4. grid position
        # 5. set
        # calculates Stokes parameters and PA/PD
        for date in dates:
            sdf = source[source['date'] == date]
            objects = sdf['d_object'].unique()
            for obj in objects:
                d_pos = sdf[sdf['d_object'] == obj]['d_grid_pos'].unique()
                sets = sdf[sdf['d_object'] == obj]['d_set'].unique()
                stars = sdf['letters'].unique()
                for star in stars:
                    obj_sdf = sdf[sdf['d_object'] == obj]
                    ap_sdf = obj_sdf[obj_sdf['aperture_set'] == aperture]
                    star_sdf = ap_sdf[ap_sdf['letters'] == star]
                    for pos in d_pos:
                        pos_sdf = star_sdf[star_sdf['d_grid_pos'] == pos]
                        for i in sets:
                            pol_sdf = pos_sdf[pos_sdf['d_set'] == i]
                            d_pols = pol_sdf['d_pol'].unique()
                            time = pol_sdf['time'].iloc[0]
                            # pol_sdf[""]
                            for pol in d_pols:
                                try:
                                    mag = pol_sdf[pol_sdf['d_pol'] == pol]['mag'].iloc[0]
                                except:
                                    print(i, pol, pos, "skipped")
                                    continue
                                try:
                                    mag_err = float(pol_sdf[pol_sdf['d_pol'] == pol]['mag_error'])
                                except:
                                    raise Exception("%d duplicates in: %s %s_%s-%s-p%s" % \
                                     (len(pol_sdf[pol_sdf['d_pol'] == pol]['mag_error']), date, obj, pos, i, pol))
                                I[pol] = np.power(10., -0.4 * mag)
                                I_err[pol] = abs(-0.921034 * np.exp(-0.921034 * mag)) * mag_err

                            for origin in [0, 180]:
                                angles = np.array((0, 45, 90, 135)) + origin
                                if not np.prod(np.in1d(angles, d_pols)):
                                    continue
                                PD, PDerr, PA, PAerr, Q, Qerr, U, Uerr = self.getPol(I[angles[0]], I[angles[1]], I[angles[2]], I[angles[3]],
                                                                                I_err[angles[0]], I_err[angles[1]], I_err[angles[2]], I_err[angles[3]])
                                row ={cols[0]: date,
                                      cols[1]: time,
                                      cols[2]: np.nan,
                                      cols[3]: obj,
                                      cols[4]: np.nan,
                                      cols[5]: np.nan,
                                      cols[6]: pos,
                                      cols[7]: np.nan,
                                      cols[8]: np.nan,
                                      cols[9]: i,
                                      cols[10]: star,
                                      cols[11]: angles[0],
                                      cols[12]: Q,
                                      cols[13]: Qerr,
                                      cols[14]: U,
                                      cols[15]: Uerr,
                                      cols[16]: PD,
                                      cols[17]: PDerr,
                                      cols[18]: PA,
                                      cols[19]: PAerr}
                                if(len(row) != len(cols)):
                                    print("WARNING: column mismatch in calcPol")
                                df_pol = df_pol.append(row, ignore_index=True)
        return df_pol

    # drops existing table to avoid overlaps
    def dropTable(self, tablename = ""):
        if tablename == "":
            raise Exception("dropTable: no table specified")
            exit()

        mydb = mysql.connector.connect(
            host="localhost",
            user="root",
            password="admin",
            database="mydatabase"
        )

        mycursor = mydb.cursor()
        cmdline = "DROP TABLE %s;" % tablename
        try:
            mycursor.execute(cmdline)
        except:
            print("no such table: %s. nothing dropped" % tablename)
        for x in mycursor:
            print(x)

        mydb.commit()
        mydb.close()

    def createTable(self, tablename = ""):
        if tablename == "":
            raise Exception("createTable: no table specified")
            exit()
        mydb = mysql.connector.connect(
            host="localhost",
            user="root",
            password="admin",
            database="mydatabase"
        )
        mycursor = mydb.cursor()
        if tablename == "poldata":
            cmdline = "CREATE TABLE %s(" % tablename +\
                        "obs INT PRIMARY KEY AUTO_INCREMENT," +\
                        "obs_date DATE," +\
                        "obs_time VARCHAR(10)," +\
                        "MJD VARCHAR(10)," +\
                        "src VARCHAR(10)," +\
                        "src_RA VARCHAR(10)," +\
                        "src_DEC VARCHAR(10)," +\
                        "position VARCHAR(2)," +\
                        "ccd_x VARCHAR(10)," +\
                        "ccd_y VARCHAR(10)," +\
                        "dset INT," +\
                        "star VARCHAR(10)," +\
                        "origin INT," +\
                        "Q FLOAT," +\
                        "Qerr FLOAT," +\
                        "U FLOAT," +\
                        "Uerr FLOAT," +\
                        "PD FLOAT," +\
                        "PDerr FLOAT," +\
                        "PA FLOAT," +\
                        "PAerr FLOAT);"

        elif tablename == "diffphot":
            cmdline = "CREATE TABLE %s(" % tablename +\
                        "obs INT PRIMARY KEY AUTO_INCREMENT," +\
                        "obs_date DATE," +\
                        "time VARCHAR(10)," +\
                        "exptime FLOAT," +\
                        "fname VARCHAR(30)," +\
                        "aperture_set VARCHAR(20)," +\
                        "target VARCHAR(10)," +\
                        "sky FLOAT," +\
                        "fwhm FLOAT," +\
                        "ellipticity FLOAT," +\
                        "value2 FLOAT," +\
                        "value3 FLOAT," +\
                        "mag FLOAT," +\
                        "mag_error FLOAT," +\
                        "status VARCHAR(30)," +\
                        "d_object VARCHAR(10)," +\
                        "d_grid_pos VARCHAR(3)," +\
                        "d_set INT," +\
                        "d_pol INT);"

        else:
            raise Exception("Unknown table type: %s. Exiting." % tablename)
            exit()
        mycursor.execute(cmdline)
        for x in mycursor:
            print(x)
        mydb.commit()
        mydb.close()

    def writeDatabase(self, df, tablename = ""):
        if tablename == "":
            raise Exception("writeDatabase: no table name specified")
            exit()
        # drop an existing table first
        self.dropTable(tablename)
        # re-create a table
        self.createTable(tablename)
        # insert into the table
        mydb = mysql.connector.connect(
            host="localhost",
            user="root",
            password="admin",
            database="mydatabase"
        )

        mycursor = mydb.cursor()
        line_counter = 0
        for i in range(len(df)):
            line_counter += 1
            cmdline = "INSERT INTO %s(" % tablename
            if tablename == "poldata":
                sqlcolnames = ["obs_date", "obs_time", "MJD", "src", "src_RA", "src_DEC", "position",
                               "ccd_x", "ccd_y", "dset", "star", \
                               "origin", "Q", "Qerr", "U", "Uerr", "PD", \
                               "PDerr", "PA", "PAerr"]
            elif tablename == "diffphot":
                sqlcolnames = ["obs_date", "time", "exptime", "fname", "aperture_set", "target", "sky", "fwhm", "ellipticity",
                               "value2", "value3", "mag", "mag_error", "status", "d_object", "d_grid_pos",
                               "d_set", "d_pol"]
            else:
                raise Exception("Unknown table type. Exiting.")
                exit()
            colnames = np.asarray(df.columns)
            for col in sqlcolnames:
                if col != sqlcolnames[-1]:
                    cmdline += col + ', '
                else:
                    cmdline += col
            cmdline += ")"
            cmdline += " VALUES ("
            special_cols = []
            if tablename == "poldata":
                special_cols = ["date", "time", 'DEC', "RA", "x", "y", "MJD", "source", "position", "star"]
            elif tablename == "diffphot":
                special_cols = ["date", "time", "name", "target", "aperture_set", "letters", "status", "d_object", "d_grid_pos"]
            else:
                raise Exception("writeDatabase: unknown table")
                exit()
            for col in colnames:
                if col in special_cols:
                    cmdline += "\'" + str(df[col][i])+ "\'" + ", "
                elif col != colnames[-1]:
                    cmdline += str(df[col][i]) + ", "
                else:
                    cmdline += str(df[col][i])
            cmdline += ");"
            #print(cmdline)
            mycursor.execute(cmdline)
        mydb.commit()
        mydb.close()

        print("write successful: %d lines" % line_counter)

    # parses a database table into a dataframe
    def readDatabase(self, tablename = ""):
        if tablename == "":
            raise Exception("readDatabase: table name not specified")
            exit()
        colnames = []
        if tablename == "poldata":
            colnames = self.getColNames("poldata")
        elif tablename == "diffphot":
            colnames = self.getColNames("diffphot")
        df_out = pd.DataFrame(columns = colnames)
        mydb = mysql.connector.connect(
            host="localhost",
            user="root",
            password="admin",
            database="mydatabase"
        )
        mycursor = mydb.cursor()
        mycursor.execute("SELECT * FROM %s;" %tablename)
        data = mycursor.fetchall()
        rows = [] # make a list of dictionaries to combine them in a dataframe
        if tablename == "poldata":
            cols = self.getColNames("poldata")
            for i in range(len(data)):
                row = {}
                for k in range(len(cols)): # read data to corresponding keys
                    row[cols[k]] = data[i][k+1]
                rows.append(row) # append dictionary list
                df_out = df_out.append(row, ignore_index=True)
        elif tablename == "diffphot":
            cols = self.getColNames("diffphot")
            for i in range(len(data)):
                row = {}
                for k in range(len(cols)):
                    row[cols[k]] = data[i][k+1]
                rows.append(row)
                df_out = df_out.append(row, ignore_index=True)
        mydb.close()

        return df_out

if __name__ == "__main__":
    # writing test 0
    db = DatabaseHandler("./data/db_test/", "diffphot")
    db.updatePolTable()
    print(db.readDatabase("poldata"))
    # db.calcPol()
    # db.fillTable("poldata")
    # print(db.getDataframe("poldata"))

    # writing test 1
    db = DatabaseHandler()
    df = db.readFolder("./data/db_test/")
    print(df)
    db.writeDatabase(df, "diffphot")
    db.updatePolTable()
    print(db.readDatabase("poldata"))
