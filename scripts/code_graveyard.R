# code graveyard

# code to help bind multiple files that have started at the same time! ####

# grab serial number from the info file
d_serial <- filter(d_info,  str_detect(param, 'Serial number')) %>%
  select(file, serial_no = value)

# bind with od
d_od <- left_join(d_od, d_serial)

# split the file name of d_od to its useful parts
d_od <- d_od %>%
  # split file name based on each _
  separate(file, c('id', 'run', 'date', 'temp'), sep = '_') %>%
  mutate(temp = parse_number(temp)) %>%
  # calculate earliest time for each file
  group_by(run, temp, date) %>%
  mutate(min_time = min(time)) %>% 
  ungroup() %>%
  # only keep certain columns
  select(run, serial_no, date, min_time, measurement_time_s, temp, wavelength_nm, well, raw_absorbance)

# work out which file for each temperature is the first one
d_times <- d_od %>%
  # keep unique instances of group columns
  select(run, serial_no, date, temp, min_time) %>%
  distinct() %>%
  # turn date column into a time object
  mutate(date2 = ymd(date)) %>%
  group_by(serial_no) %>%
  arrange(date2) %>%
  # create order of file creation
  mutate(order = 1:n()) %>%
  ungroup() %>%
  arrange(serial_no, date)
# serial number 1550-802208 timestamp is messed up, does not even change between runs so need to check that

# calculate max measurement time for each file
d_max_time <- d_od %>%
  group_by(run, serial_no, date, temp) %>%
  summarise(max_time = max(measurement_time_s), .groups = 'drop')

# merge together
d_times <- left_join(d_times, d_max_time) %>%
  group_by(serial_no) %>%
  # calculate time to add to follow on file
  mutate(to_add = cumsum(lag(max_time, n =1, default = 0))) %>%
  # add on another five minutes as it starts at 0
  mutate(to_add = to_add + (300*(order - 1))) 

# keep only columns that are important
d_times <- select(d_times, serial_no, date, temp, to_add)

# merge d_times with d_od 
d_od <- left_join(d_od, d_times) %>%
  # add time to make sure measurement times are correct
  mutate(measurement_time_s = measurement_time_s + to_add) 

head(d_od)

# CODE GRAVEYARD ####

# calculate measurement time in seconds from the start time - this does not seem to work as the internal clocks of the readers seem to be shit!
# if all times were correct on the plate readers, this would be the optimal way, but it does not seem to be the case, even relatively
#d_od <- group_by(d_od, run, well, serial_no) %>%
#mutate(measurement_time_s = time - min(time),
#       measurement_time_s = as.numeric(measurement_time_s)) %>%
#arrange(measurement_time_s) %>%
#ungroup()
