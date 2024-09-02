% Save deformed configuration
outFile = 'Output/xnCompressed.out';
if jsave == 1      
    if jcont == 1  
  
        fid = fopen(outFile, 'w');
        for js = 1:nCP
            fprintf(fid, '%18.9e %18.9e \n', xn(js, :));
        end
        fclose(fid);        
        
    end
end 