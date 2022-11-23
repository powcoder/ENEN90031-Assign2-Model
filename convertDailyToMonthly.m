function data_monthly = convertDailyToMonthly( data )
%convertDailyToMonthly agregates daily data to monthly data
%   
%   Written by Tim Peterson, University of Melbourne, May 2010


    nrows = size(data,1);            
    month_row = 1;    
    currentMonth = data( 1, 2);
    currentYear = data(1,1);
    data_monthly = [ data(1,1), data(1,2), eomday( data(1,1), data(1,2)), data(1,4:end)];
    daysAdded = 1;
    
    for day_row=2: nrows        
       
       if data(day_row,2) ==  currentMonth;
           data_monthly( month_row, 4:end) = data_monthly( month_row, 4:end) + data( day_row , 4:end);
           daysAdded=daysAdded+1;
       else
           %check added by Andrew Western
           if daysAdded ~= eomday(currentYear,currentMonth)
               data_monthly( month_row, 4:end) = nan(size(data_monthly( month_row, 4:end)));
           end
           month_row = month_row+1;
           currentMonth = data( day_row, 2);
           currentYear = data(day_row,1);
           data_monthly(month_row,:) = [ data(day_row,1), data(day_row,2), eomday( data(day_row,1), data(day_row,2)), data(day_row,4:end)];
           daysAdded = 1;
       end        
        
    end
    
    if daysAdded ~= eomday(currentYear,currentMonth)
        data_monthly= data_monthly(1:end-1, :);
    end
    
    if any(isnan(data_monthly(1, 4:end)));
        data_monthly= data_monthly(2:end, :);
    end

end

