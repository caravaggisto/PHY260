import datetime
import os

# this could be set in a file for ALL homeworks...
local_data_directory = 'data'
local_figs_directory = 'figs'

# this is set for particular homeworks
figs_directory = os.getcwd() + '/' + local_figs_directory
data_directory = os.getcwd() + '/' + local_data_directory

def get_newest_file_by_type_in_dir( file_extension_including_period, local_data_dir = local_data_directory ):
    # taken from http://ubuntuforums.org/showthread.php?t=1526010
    os.chdir( local_data_dir )
    # change into the data directory
    filelist = [filename for filename in os.listdir( os.getcwd() ) if filename.endswith( file_extension_including_period ) ]
    # make a list of all files and folders in the current working directory
    filelist = filter( lambda x: not os.path.isdir( x ), filelist )
    # get rid of folders in the list. we only want files.
    newest = max( filelist, key = lambda x: os.stat( x ).st_mtime )
    # get newest file in the list
    path_to_newest = os.getcwd() + '/' + newest
    os.chdir( '..' )
    # change back to the parent directory.
    return path_to_newest

def filepath_for_now( local_folder_for_file = 'data', hw_num = 6 ):
    ''' args:
        local_folder_for_file = 'data' or 'figs'
        hw_num = N for hw #N, e.g. hw #6 --> hw_num = 6.
    NOTE: this return the path relative to the working directory of the *calling* function!
    '''
    # return something like /path/to/workspace/projectname/hw6/data/2012-04-02-15:44:02.p
    time_string = datetime.datetime.now().strftime( "%Y-%m-%d-%H:%M:%S" )
    # taken from http://stackoverflow.com/questions/6327498/
    full_path_to_folder_to_file = os.getcwd() + '/' + local_folder_for_file
    if os.path.exists( full_path_to_folder_to_file ):
        return full_path_to_folder_to_file + '/' + time_string
    else:
        raise UserWarning( "the path \n\t" + full_path_to_folder_to_file + "\ndoesn't exist!!" )
