 function Distance=DistanceToBorder(Position,Corner_position)

%     Position=FullNeuronBehaviorDataSet.FullNeuronBehaviorMatrix{1,1}(:,2:3);
%     Corner_position=FullNeuronBehaviorDataSet.EBC_results{1,4}.QP;
%%
     for i=1:1:size(Position,1)

         Position_now=[squeeze(Position(i,:)),0];
         V1=[Corner_position(1,:),0];
         V2=[Corner_position(2,:),0];
         V3=[Corner_position(3,:),0];
         V4=[Corner_position(4,:),0];

        a = V1 - V2;
        b = Position_now - V2;
        d1 = norm(cross(a,b)) /norm(a);


        a = V2 - V3;
        b = Position_now - V3;
        d2 = norm(cross(a,b)) /norm(a);


        a = V3 - V4;
        b = Position_now - V4;
        d3 = norm(cross(a,b)) /norm(a);


        a = V4 - V1;
        b = Position_now - V1;
        d4 = norm(cross(a,b)) /norm(a);

        Distance(i)=min([d1,d2,d3,d4]);
    end


 end