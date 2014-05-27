function [  ] = interactions_to_table( interactions_all, time_range, filename )
%Take interactions from MATLAB form (interactions_all) and write them to an
%R-readable table form.

t_start = time_range(1);
t_end = time_range(2);

N = size(interactions_all,1);

fid=fopen(filename,'w');

fprintf(fid,'day,leader,follower,event.type,t1,t2,t3,strength,disp,lead.dx.12,lead.dy.12,foll.dx.12,foll.dy.12,lead.dx.23,lead.dy.23,foll.dx.23,foll.dy.23,dyad.dx.t2,dyad.dy.t2\n')
for i = 1:N
    i
    for j = (i+1):N
        interactions = interactions_all{i,j};
        for k = 1:length(interactions)
            if interactions(k).time_idxs(1) >= t_start && interactions(k).time_idxs(3) <= t_end
                fprintf(fid,[num2str(interactions(k).day) ',' num2str(interactions(k).leader) ',' num2str(interactions(k).follower) ',' ...
                    interactions(k).type ',' num2str(interactions(k).time_idxs(1)) ',' num2str(interactions(k).time_idxs(2)) ',' num2str(interactions(k).time_idxs(3)) ',' ...
                    num2str(interactions(k).strength_mult) ',' num2str(interactions(k).disp_mult) ',' ...
                    num2str(interactions(k).leader_disp_12(1)) ',' num2str(interactions(k).leader_disp_12(2)) ',' ...
                    num2str(interactions(k).follower_disp_12(1)) ',' num2str(interactions(k).follower_disp_12(2)) ',' ...
                    num2str(interactions(k).leader_disp_23(1)) ',' num2str(interactions(k).leader_disp_23(2)) ',' ...
                    num2str(interactions(k).follower_disp_23(1)) ',' num2str(interactions(k).follower_disp_23(2)) ...
                    num2str(interactions(k).dyad_vector_t2(1)) ',' num2str(interactions(k).dyad_vector_t2(2)) ...
                    '\n']);
            end
        end
    end
end

fclose(fid)


end

