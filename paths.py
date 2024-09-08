import os

# List of directory and files paths
current_dir = os.path.dirname(os.path.abspath(__file__))

#data
data_path = os.path.join(current_dir, 'data')
#   Juno
juno_path = os.path.join(data_path, 'Juno')
#   Izumo
izumo_path = os.path.join(data_path, 'Izumo')


#conservation_analysis
conservation_analysis_path = os.path.join(current_dir, 'conservation_analysis')
#   Juno
juno_conservation_path = os.path.join(conservation_analysis_path, 'Juno')
#   Izumo
izumo_conservation_path = os.path.join(conservation_analysis_path, 'Izumo')
