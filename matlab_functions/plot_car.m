function [car,V]=plot_car()
    car1=stlread('car_1.stl');
    V=car1.Points/5;
    F=car1.ConnectivityList;
    V(:,3)=-V(:,3);
    V(:,1)=-V(:,1);
    car=patch('Vertices',V,'Faces',F);
    car.FaceColor='r';
    car.FaceAlpha=0.1;
    car.EdgeAlpha=0;
    car.PickableParts='none';
end