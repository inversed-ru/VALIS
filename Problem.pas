{ 
Copyright (c) Peter Karpov 2010 - 2017.

Usage of the works is permitted provided that this instrument is retained with 
the works, so that any entity that uses the works is notified of this instrument.

DISCLAIMER: THE WORKS ARE WITHOUT WARRANTY.
}
{$IFDEF FPC} {$MODE DELPHI} {$ENDIF}
unit Problem; ///////////////////////////////////////////////////////////////////////
{
>> Description
   VALIS problem definition unit for real-valued classification problems.

>> Author
   Peter Karpov
   Email    : PeterKarpov@inversed.ru
   Homepage : inversed.ru
   GitHub   : inversed-ru
   Twitter  : @inversed_ru
   
>> ToDo
   - Examine various crossover and mutation operators employed in real coded GAs
   - Initialization option: apply crossover to the training data
   ? Covariance matrix adaptation
   
>> Changelog
   1.4  : 2017.12.16 ~ InvLibs compatibility
                     ~ Cleanup
   1.3  : 2017.02.09 + Intermediate and midpoint crossovers
                     ~ Modified initialization to exclude duplicates
                     * AssignAntibody not copying some attributes
                     + Various radius initialization methods
   1.2  : 2017.01.03 * Incorrect metric inheritance during crossover
   1.1  : 2016.12.16 * Incorrect mutation
   1.0  : 2016.12.15 ~ Rewritten from scratch
   0.0  : 2009.05.21 + Initial version
   Notation: + added, - removed, * fixed, ~ changed

}
interface ///////////////////////////////////////////////////////////////////////////
uses
      Arrays;
type
      TCrossoverType = (ctBinary, ctIntermediate, ctMidpoint);
      TRadiusInit    = 
      (  riMean, riMeanUniform, riMeanTriangle, riMeanExp, 
         riDistToNearest, riDistToRandom
      );
const
      AdaptRadius       =  True;
      AdaptDimWeights   =  False;
      AdaptMetric       =  False;
      RandInit          =  False;
      RadiusInit        =  riDistToRandom;
      CrossoverType     =  ctBinary;

type
      TAntigen =
         record
         x        :  TRealArray;
         Class_   :  Integer;
         end;

      TAntigens =
         record
         _        :  array of TAntiGen;
         N, NVars,
         NClasses :  Integer;
         end;

      TClassDistr = TRealArray;

      TAntibody =
         record
         // Generic fields
         SumW,
         SumWCorrect,
         Accuracy,
         KShare,
         Fitness        :  Real;
         MajorityClass  :  Integer;
         ClassDistr     :  TRealArray;

         // Problem-specific fields
         // #TODO should be a separate type?
         x, DimWeight   :  TRealArray;
         Radius,
         Metric         :  Real;
         end;

      TAntibodies =
         record
         _              :  array of TAntibody;
         OrderFitness,
         OrderSumW      :  TIntArray;
         Accuracy,
         Coverage       :  Real;
         N              :  Integer;
         end;
         
{-----------------------<< Antigens >>----------------------------------------------}

// Create a new empty array of Antigens
procedure InitAntigens(
   var   Antigens    :  TAntigens);
         overload;

// Create a new empty array of antigens, set the number of variables and classes
procedure InitAntigens(
   var   Antigens    :  TAntigens;
         NVars,
         NClasses    :  Integer);
         overload; 

// Append an Antigen to the Antigens
procedure AppendAntigen(
   var   Antigens :  TAntigens;
   const Antigen  :  TAntigen); 
   
// Load the Data from Path
procedure LoadData(
   var   Data     :  TAntigens;
   const Path     :  AnsiString);
   
{-----------------------<< Antibodies >>--------------------------------------------}

// ABTo := ABFrom
procedure AssignAntibody(
   var   ABTo     :  TAntibody;
   const ABFrom   :  TAntibody);

// ABsTo := ABsFrom
procedure CopyAntibodies(
   var   ABsTo    :  TAntibodies;
   const ABsFrom  :  TAntibodies);           

// Return the distance between the Antigen and the Antibody
function Distance(
   const Antigen  :  TAntigen;
   const Antibody :  TAntibody
         )        :  Real;  

// Create a Child by crossing over ParentA and ParentB
procedure Crossover(
   var   Child    :  TAntibody;
   const ParentA,
         ParentB  :  TAntibody);

// Mutate the Antibody using a given mutation Rate
procedure Mutate(
   var   Antibody    :  TAntibody;
         Rate        :  Real);

// Return the Order of Antigens sorted by distance to x
// #HACK should not be in the problem unit
procedure SortByDistance(
   var   Order    :  TIntArray;
   const Antigens :  TAntigens;
   const x        :  TRealArray);         

// Create a population of N Antibodies based on Antigens
procedure InitAntibodies(
   var   Antibodies  :  TAntibodies;
         N           :  Integer;
   const Antigens    :  TAntigens);  

// Save antibody parameters to Path
procedure SaveAntibodiesStats(
   const Path           :  AnsiString;
   const Antibodies     :  TAntibodies);

// Save antibody parameters required for visualization to Path
procedure SaveAntibodiesVis(
   const Path           :  AnsiString;
   const Antibodies     :  TAntibodies;
   const Antigens       :  TAntigens);      

implementation //////////////////////////////////////////////////////////////////////
uses
      InvSys,
      Math,
      ExtraMath,
      Statistics,
      StringUtils,
      RandVars,
      Sorting;

{-----------------------<< Antigens >>----------------------------------------------}

// AGTo := AGFrom
procedure Assign(
   var   AGTo     :  TAntigen;
   const AGFrom   :  TAntigen);
   begin
   AGTo.Class_ := AGFrom.Class_;
   AGTo.x := Copy(AGFrom.x);
   end;


// Create a new empty array of Antigens
procedure InitAntigens(
   var   Antigens    :  TAntigens);
         overload;
   begin
   SetLength(Antigens._, 0);
   Antigens.N := 0;
   end;


// Create a new empty array of antigens, set the number of variables and classes
procedure InitAntigens(
   var   Antigens    :  TAntigens;
         NVars,
         NClasses    :  Integer);
         overload;
   begin
   InitAntigens(Antigens);
   Antigens.NVars := NVars;
   Antigens.NClasses := NClasses;
   end;


// Append an Antigen to the Antigens
procedure AppendAntigen(
   var   Antigens :  TAntigens;
   const Antigen  :  TAntigen);
   begin
   with Antigens do
      begin
      Inc(N);
      if N > Length(Antigens._) then
         SetLength(Antigens._, 2 * Antigens.N);
      Assign(Antigens._[N - 1], Antigen);
      end;
   end;


// Load the Data from Path
procedure LoadData(
   var   Data     :  TAntigens;
   const Path     :  AnsiString);
   var
         i        :  Integer;
         S        :  AnsiString;
         Antigen  :  TAntigen;
         FileData :  Text;
   begin
   with Data do
      begin
      // Read the header
      OpenRead(FileData, Path);
      ReadLn(FileData, S);
      NVars := SubstrCount(S, Tab);

      // Read the data
      InitAntigens(Data);
      NClasses := 1;
      SetLength(Antigen.x, NVars);
      repeat
         for i := 0 to NVars - 1 do
            Read(FileData, Antigen.x[i]);
         ReadLn(FileData, Antigen.Class_);
         if Antigen.Class_ >= NClasses then
            Inc(NClasses);
         AppendAntigen(Data, Antigen);
      until EoF(FileData);
      Close(FileData);
      end;
   end;

{-----------------------<< Antibodies >>--------------------------------------------}
   
// ABTo := ABFrom
procedure AssignAntibody(
   var   ABTo     :  TAntibody;
   const ABFrom   :  TAntibody);
   begin
   // Problem-specific fields
   ABTo.x               := Copy(ABFrom.x);
   ABTo.DimWeight       := Copy(ABFrom.DimWeight);
   ABTo.Radius          :=      ABFrom.Radius;
   ABTo.Metric          :=      ABFrom.Metric;

   // Generic fields
   ABTo.SumW            :=      ABFrom.SumW;
   ABTo.SumWCorrect     :=      ABFrom.SumWCorrect;
   ABTo.Accuracy        :=      ABFrom.Accuracy;
   ABTo.KShare          :=      ABFrom.KShare;
   ABTo.Fitness         :=      ABFrom.Fitness;
   ABTo.MajorityClass   :=      ABFrom.MajorityClass;
   ABTo.ClassDistr      := Copy(ABFrom.ClassDistr);
   end;


// ABsTo := ABsFrom
procedure CopyAntibodies(
   var   ABsTo    :  TAntibodies;
   const ABsFrom  :  TAntibodies);
   var
         i        :  Integer;
   begin
   ABsTo.N            :=      ABsFrom.N;
   ABsTo.Accuracy     :=      ABsFrom.Accuracy;
   ABsTo.Coverage     :=      ABsFrom.Coverage;
   ABsTo.OrderFitness := Copy(ABsFrom.OrderFitness);
   ABsTo.OrderSumW    := Copy(ABsFrom.OrderSumW);
   SetLength(ABsTo._, ABsTo.N);
   for i := 0 to ABsTo.N - 1 do
      AssignAntibody(ABsTo._[i], ABsFrom._[i]);
   end;
   
   
// Root mean square based measure of distance between A and B
function DistanceRMS(
   const A, B  :  TRealArray
   )           :  Real;
   var
         i     :  Integer;
         Sum   :  Real;
   begin
   Sum := 0;
   for i := 0 to Length(A) - 1 do
      Sum := Sum + Sqr(A[i] - B[i]);
   Result := Sqrt(Sum / Length(A));
   end;
   

// Return the distance between the Antigen and the Antibody
function Distance(
   const Antigen  :  TAntigen;
   const Antibody :  TAntibody
         )        :  Real;
   var
         i, NVars :  Integer;
         d        :  TRealArray;
   begin
   NVars := Length(Antibody.x);
   SetLength(d, NVars);
   for i := 0 to NVars - 1 do
      d[i] := Power(Abs(Antigen.x[i] - Antibody.x[i]), Antibody.Metric);
   Result := Power(WeightedAverage(d, Antibody.DimWeight), 1 / Antibody.Metric);
   end;


// Return a random blending factor depending on CrossoverType
function RandCrossAlpha :  Real;
   begin
   case CrossoverType of
      ctBinary       :  Result := Random(2);
      ctIntermediate :  Result := Random;
      ctMidpoint     :  Result := 1 / 2;
      else              Result := 1 / 2;
                        Assert(False);
      end;
   end;
   

// Create a Child by crossing over ParentA and ParentB
// #HACK need a parameter iterator?
procedure Crossover(
   var   Child    :  TAntibody;
   const ParentA,
         ParentB  :  TAntibody);
   var
         i, NVars :  Integer;
   begin
   // Center coordinates
   NVars := Length(ParentA.x);
   SetLength(Child.x, NVars);
   for i := 0 to NVars - 1 do
      Child.x[i] := Blend(ParentA.x[i], ParentB.x[i], RandCrossAlpha);

   // Dimension weights
   SetLength(Child.DimWeight, NVars);
   for i := 0 to NVars - 1 do
      Child.DimWeight[i] := Blend(ParentA.DimWeight[i], ParentB.DimWeight[i], RandCrossAlpha);
   Normalize(Child.DimWeight, normUnitSum);

   // Radius and metric
   Child.Radius := Blend(ParentA.Radius, ParentB.Radius, RandCrossAlpha);
   Child.Metric := Blend(ParentA.Metric, ParentB.Metric, RandCrossAlpha);
   end;


// Mutate the Antibody using a given mutation Rate
// #HACK need a parameter iterator?
procedure Mutate(
   var   Antibody    :  TAntibody;
         Rate        :  Real);
   var
         i, NVars    :  Integer;
         Changed,
         DimWChanged :  Boolean;
   const
         MinStep     =  1 / 256;
         MaxStep     =  1;
         RScale      =  1 / 4;
   begin
   Changed := False;
   DimWChanged := False;
   NVars := Length(Antibody.x);
   repeat
      // Center coordinates
      for i := 0 to NVars - 1 do
         if Random < Rate then
            begin
            Changed := True;
            Antibody.x[i] := Antibody.x[i] + RandSign * LogBlend(MinStep, MaxStep, Random);
            end;

      // Dimension weights
      if AdaptDimWeights then
         for i := 0 to NVars - 1 do
            if Random < Rate then
               begin
               Changed := True;
               DimWChanged := True;
               Antibody.DimWeight[i] := Antibody.DimWeight[i] * Exp(RScale * RandGauss);
               end;
      if DimWChanged then
         Normalize(Antibody.DimWeight, normUnitSum);

      // Radius
      if (Random < Rate) and AdaptRadius then
         begin
         Changed := True;
         Antibody.Radius := Antibody.Radius * Exp(RScale * RandGauss);
         end;

      // Metric
      if (Random < Rate) and AdaptMetric then
         begin
         Changed := True;
         Antibody.Metric := Exp(Sqrt(Ln(2) / 2) * RandGauss + Ln(2));
         end;
   until Changed;
   end;


// Return the Order of Antigens sorted by L2 distance to x
// #HACK should not be in the problem unit
procedure SortByDistance(
   var   Order    :  TIntArray;
   const Antigens :  TAntigens;
   const x        :  TRealArray);
   var
         i        :  Integer;
         D        :  TRealArray;
   begin
   SetLength(D, Antigens.N);
   for i := 0 to Antigens.N - 1 do
      D[i] := DistanceRMS(Antigens._[i].x, x);
   OrderRealArray(Order, D, soAScending);
   end;


// Create a population of N Antibodies based on Antigens
procedure InitAntibodies(
   var   Antibodies  :  TAntibodies;
         N           :  Integer;
   const Antigens    :  TAntigens);
   var
         i, j, 
         NVars       :  Integer;
         Order,
         DistOrder   :  TIntArray;
         MeanDist    :  Real;
   begin
   Antibodies.N := N;
   SetLength(Antibodies._, N);
   NVars := Antigens.NVars;
   RandPerm(Order, Antigens.N, {Base:} 0);
   
   // Calculate the mean pairwise antigen distance if necessary
   if RadiusInit in [riMean, riMeanUniform, riMeanTriangle, riMeanExp] then
      begin
      MeanDist := 0;
      for i := 0 to Antigens.N - 1 do
         for j := 0 to i - 1 do
            MeanDist := MeanDist + DistanceRMS(Antigens._[i].x, Antigens._[j].x);
      MeanDist := MeanDist / NPairs(Antigens.N);
      end;

   // Initialize the antobodies
   for i := 0 to N - 1 do
      with Antibodies._[i] do
         begin
         if RandInit then
            begin
            SetLength(x, NVars);
            for j := 0 to NVars - 1 do
               x[j] := RandGauss;
            end
         else
            x := Copy(Antigens._[Order[i mod Antigens.N]].x);
         case RadiusInit of
            riMean:
               Radius := MeanDist;
            riMeanUniform:
               Radius := MeanDist * Random * 2;
            riMeanTriangle:
               Radius := MeanDist * (Random + Random);
            riMeanExp:
               Radius := MeanDist * RandExp;
            riDistToNearest:
               begin
               SortByDistance(DistOrder, Antigens, x);
               j := 1;
               Radius := 0;
               while Radius = 0 do
                  begin
                  Radius := DistanceRMS(x, Antigens._[DistOrder[j]].x);
                  Inc(j);
                  end;
               end;
            riDistToRandom:
               begin
               Radius := 0;
               while Radius = 0 do
                  Radius := DistanceRMS(x, Antigens._[Random(Antigens.N)].x);
               end;
            else
               Radius := 0;
               Assert(False);
            end;

         // #TODO Randomize DimWeights and Metric if adaptive?
         Metric := 2;
         SetLength(DimWeight, NVars);
         for j := 0 to NVars - 1 do
            DimWeight[j] := 1;
         end;
   end;
   

// Save antibody parameters to Path
procedure SaveAntibodiesStats(
   const Path           :  AnsiString;
   const Antibodies     :  TAntibodies);
   var
         i, j           :  Integer;
         FileAntibodies :  Text;
   begin
   OpenWrite(FileAntibodies, Path);
   WriteLn(FileAntibodies,
      'Id'        ,  Tab,
      'Fitness'   ,  Tab,
    //'Accuracy'  ,  Tab,
    //'NShare '   ,  Tab,
      'NBound'    ,  Tab,
      'Radius'    ,  Tab,
      'Shape'     );
   for i := 0 to Antibodies.N - 1 do
      begin
      j := Antibodies.OrderFitness[Antibodies.N - 1 - i];
      WriteLn(FileAntibodies,
         i                                      ,  Tab,
         Antibodies._[j].Fitness                ,  Tab,
       //Antibodies._[j].Accuracy               ,  Tab,
       //1 / Antibodies._[j].KShare             ,  Tab,
         Antibodies._[j].SumW                   ,  Tab,
         Antibodies._[j].Radius                 ,  Tab,
         Antibodies._[j].Metric                 );
      end;
   Close(FileAntibodies);
   end;
   
   
// Save antibody parameters required for visualization to Path
procedure SaveAntibodiesVis(
   const Path           :  AnsiString;
   const Antibodies     :  TAntibodies;
   const Antigens       :  TAntigens);
   var
         i, j           :  Integer;
         FileAntibodies :  Text;
   begin
   OpenWrite(FileAntibodies, Path);
   for i := 0 to Antibodies.N - 1 do
      begin
      // Coordinates
      for j := 0 to High(Antibodies._[i].x) do
         Write(FileAntibodies, Antibodies._[i].x[j], Tab);

      // Class distribution
      for j := 0 to High(Antibodies._[i].ClassDistr) do
         Write(FileAntibodies, Antibodies._[i].ClassDistr[j], Tab);

      // Radius
      Write(FileAntibodies, Antibodies._[i].Radius);

      // Indices of bound antigens
      for j := 0 to Antigens.N - 1 do
         begin
         Write(FileAntibodies, Tab);
         if (Distance(Antigens._[j], Antibodies._[i]) /
             Antibodies._[i].Radius) <= 1 then
            Write(FileAntibodies, 1 + j);
         end;
      
      WriteLn(FileAntibodies);
      end;
   Close(FileAntibodies);
   end;

end.
