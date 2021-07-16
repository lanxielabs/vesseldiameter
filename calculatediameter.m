function [vec_width1,time_vec] = calculatediameter(file_sc,numtrial,numframes,pos)



files = dir_sorted(fullfile(file_sc,'*.sc'));

%%number of trials

zlength = size(pos,3);
d1 = zeros(numtrial,numframes,zlength);







for trial = 1:numtrial
    %%number of sc in a trial is 162
    SCt = read_subimage(files,-1,-1,trial:numframes*trial);
    
    
    for frame = 1:numframes
            SC2 = SCt(:,:,frame);
            SC2 = SC2';
            %imshow(SC2,[0.02 0.4]);
            for i=1:zlength
                
            v1 = pos(:,:,i);
            
            count1 = 1;
            
            dist = pdist(v1);
            if dist(2)>dist(1)
            
                if v1(2,1)<v1(3,1)
                start1 = v1(2,1);
                finish1 = v1(3,1);
                else 
                start1 = v1(3,1);
                finish1 = v1(2,1);
                end
            
            
                for t =start1:1:finish1
                    c1 = improfile(SC2,[t t+v1(1,1)-v1(2,1)],[v1(2,2)+(t-v1(2,1))/(v1(3,1)-v1(2,1))*(v1(3,2)-v1(2,2)) v1(1,2)+(t-v1(2,1))/(v1(3,1)-v1(2,1))*(v1(3,2)-v1(2,2))]);
                    tmp1(count1) = calculatedistance(c1);
                    count1 = count1+1;
                end
            
            
           
            
                
            else
                if v1(1,1)<v1(2,1)
                    start1 = v1(2,1);
                    finish1 = v1(1,1);
                else
                    start1 = v1(1,1);
                    finish1 = v1(2,1);
                end
                for t = start1:1:finish1
                    c1 = improfile(SC2,[t t+v1(4,1)-v1(1,1)],[v1(1,2)+(t-v1(1,1))/(v1(2,1)-v1(1,1))*(v1(2,2)-v1(1,2)) v1(4,2)+(t-v1(1,1))/(v1(2,1)-v1(1,1))*(v1(2,2)-v1(1,2))]);
                    tmp1(count1) = calculatedistance(c1);
                    count1 = count1+1;
                end
               
            end
            d1(trial,frame,i)=mean(tmp1);
            end
    end
end

vec_width1 = mean(d1,1);


time_vec= 0.1 * [0:numframes-1];        %in (sec)