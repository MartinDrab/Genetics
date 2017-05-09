Unit MainForm;

Interface

Uses
  Winapi.Windows, Winapi.Messages, System.SysUtils, System.Variants,
  System.Classes, Vcl.Graphics,
  Vcl.Controls, Vcl.Forms, Vcl.Dialogs, Vcl.ExtCtrls, Vcl.StdCtrls,
  Vcl.ComCtrls, Generics.Collections,
  VCFRecord;

Type
  EVCFIncidentType = (
    vcfitTruePositive,
    vcfitFalseNegative,
    vcfitFalsePositive
  );
  EVCFIncidentTypeSet = Set Of EVCFIncidentType;

  TVCFIncident = Record
    IncidentType : EVCFIncidentType;
    VCFRecord : TVCFRecord;
    end;

  TMainFrm = class(TForm)
    MainPanel: TPanel;
    SetsAndFilesGroupBox: TGroupBox;
    SettingsGroupBox: TGroupBox;
    ResultListView: TListView;
    TruthSetFileEdit: TEdit;
    TestSetFileEdit: TEdit;
    Label1: TLabel;
    Label2: TLabel;
    TruthSetBrowseButton: TButton;
    TestSetBrowseButton: TButton;
    SNPsCheckBox: TCheckBox;
    InsertionsCheckBox: TCheckBox;
    DeletionsCheckBox: TCheckBox;
    TruePositivesCheckBox: TCheckBox;
    FalseNegativesCheckBox: TCheckBox;
    FalsePositivesCheckBox: TCheckBox;
    TruthSetOpenDialog: TOpenDialog;
    TestSetOpenDialog: TOpenDialog;
    LowerPanel: TPanel;
    LoadButton: TButton;
    HideButton: TButton;
    TPsLabel: TLabel;
    FNsLabel: TLabel;
    FPsLabel: TLabel;
    MinAWFilterEdit: TEdit;
    MinAQFilterEdit: TEdit;
    Label3: TLabel;
    Label4: TLabel;
    procedure FormCreate(Sender: TObject);
    procedure FormDestroy(Sender: TObject);
    procedure FilterCheckBoxClick(Sender: TObject);
    procedure HideButtonClick(Sender: TObject);
    procedure BrowseButtonClick(Sender: TObject);
    procedure LoadButtonClick(Sender: TObject);
    procedure ResultListViewData(Sender: TObject; Item: TListItem);
    procedure ResultListViewAdvancedCustomDrawItem(Sender: TCustomListView;
      Item: TListItem; State: TCustomDrawState; Stage: TCustomDrawStage;
      var DefaultDraw: Boolean);
  Private
    FTruthFile : TVCFFile;
    FTestFile : TVCFFile;
    FFilteredTruth : TList<TVCFRecord>;
    FFilteredTest : TList<TVCFRecord>;
    FFilterSet : EVCFRecordTypeSet;
    FFilterSetIncidents : EVCFIncidentTypeSet;
    FIncidents : TList<TVCFIncident>;

    FMinAWFilter : Cardinal;
    FMinAQFilter : Cardinal;

    FTPsCount : Cardinal;
    FFNsCount : Cardinal;
    FFPsCount : Cardinal;

    Procedure BuildFilterSet;
    Procedure ApplyFilterSet(AFile:TVCFFile; AList:TList<TVCFRecord>);
    Procedure ProcessIncidents;
    Procedure AddIncident(AType:EVCFIncidentType; ARecord:TVCFRecord; AList:TList<TVCFIncident>);
  end;

Var
  MainFrm: TMainFrm;

Implementation

Procedure TMainFrm.BuildFilterSet;
begin
FFilterSet := [vcfrtReplace];
If SNPsCheckBox.Checked Then
  FFilterSet := FFilterSet + [vcfrtSNP];

If DeletionsCheckBox.Checked Then
  FFilterSet := FFilterSet + [vcfrtDeletion];

If InsertionsCheckBox.Checked Then
  FFilterSet := FFilterSet + [vcfrtInsertion];

FFilterSetIncidents := [];
If TruePositivesCheckBox.Checked Then
  FFilterSetIncidents := FFilterSetIncidents + [vcfitTruePositive];

If FalseNegativesCheckBox.Checked Then
  FFilterSetIncidents := FFilterSetIncidents + [vcfitFalseNegative];

If FalsePositivesCheckBox.Checked Then
  FFilterSetIncidents := FFilterSetIncidents + [vcfitFalsePositive];

FMinAWFilter := StrToInt64(MinAWFilterEdit.Text);
FMinAQFilter := StrToInt64(MinAQFilterEdit.Text);
end;

Procedure TMainFrm.ApplyFilterSet(AFile:TVCFFile; AList:TList<TVCFRecord>);
begin
If Assigned(AFile) Then
  begin
  AList.Clear;
  AFile.Filter(FFilterSet, AList);
  end;
end;

Procedure TMainFrm.BrowseButtonClick(Sender: TObject);
Var
  e : TEdit;
  od : TOpenDialog;
begin
If Sender = TruthSetBrowseButton Then
  begin
  e := TruthSetFIleEdit;
  od := TruthSetOpenDialog;
  end
Else If Sender = TestSetBrowseButton Then
  begin
  e := TestSetFIleEdit;
  od := TestSetOpenDialog;
  end
Else Raise Exception.Create('Invalid brose button');

If od.Execute Then
  e.Text := od.FileName;
end;

Procedure TMainFrm.AddIncident(AType:EVCFIncidentType; ARecord:TVCFRecord; AList:TList<TVCFIncident>);
Var
  incident : TVCFIncident;
begin
If (AType In FFilterSetIncidents) Then
  begin
  Case AType Of
    vcfitTruePositive : Inc(FTPsCount);
    vcfitFalseNegative: Inc(FFNsCount);
    vcfitFalsePositive: Inc(FFPsCount);
    end;

  incident.IncidentType := AType;
  incident.VCFRecord := ARecord;
  AList.Add(incident);
  end;
end;

Procedure TMainFrm.ProcessIncidents;
Var
  I : Integer;
  incident : TVCFIncident;
  r1, r2 : TVCFRecord;
  index1, index2 : Integer;
  newList : TList<TVCFIncident>;
begin
I := 0;
WHile (I < FFilteredTest.Count) Do
  begin
  r2 := FFilteredTest[I];
  If r2.AW < FMinAWFilter Then
    FFilteredTest.Delete(I)
  Else Inc(I);
  end;

FTPsCount := 0;
FFNsCount := 0;
FFPsCount := 0;
newList := TList<TVCFIncident>.Create;
index1 := 0;
index2 := 0;
While (index1 < FFilteredTruth.Count) And (index2 < FFilteredTest.Count) Do
  begin
  r1 := FFilteredTruth[index1];
  r2 := FFilteredTest[index2];
  If TVCFRecord.Equal(r1, r2) Then
    begin
    AddIncident(vcfitTruePositive, r2, newList);
    Inc(index1);
    Inc(index2);
    end
  Else If r1.Pos < r2.Pos Then
    begin
    AddIncident(vcfitFalseNegative, r1, newList);
    Inc(index1);
    end
  Else If r1.Pos > r2.Pos Then
    begin
    AddIncident(vcfitFalsePositive, r2, newList);
    Inc(index2);
    end
  Else begin
    AddIncident(vcfitFalseNegative, r1, newList);
    AddIncident(vcfitFalsePositive, r2, newList);
    Inc(index1);
    Inc(index2);
    end;
  end;

While index1 < FFilteredTruth.Count Do
  begin
  r1 := FFilteredTruth[index1];
  AddIncident(vcfitFalseNegative, r1, newList);
  Inc(index1);
  end;

While index2 < FFilteredTest.Count Do
  begin
  r2 := FFilteredTest[index2];
  AddIncident(vcfitFalsePositive, r2, newList);
  Inc(index2);
  end;

ResultListVIew.Items.Count := 0;
FIncidents.Free;
FIncidents := newList;
ResultListView.Items.Count := FIncidents.Count;
TPsLabel.Caption := Format('TPs: %u', [FTPsCount]);
FNsLabel.Caption := Format('FNs: %u', [FFNsCount]);
FPsLabel.Caption := Format('FPs: %u', [FFPsCount]);
end;

Procedure TMainFrm.ResultListViewAdvancedCustomDrawItem(Sender: TCustomListView;
  Item: TListItem; State: TCustomDrawState; Stage: TCustomDrawStage;
  var DefaultDraw: Boolean);
Var
  c : TColor;
  i : TVCFIncident;
begin
i := FIncidents[Item.Index];
ResultListView.Canvas.Font.Style := [fsBold];
Case i.IncidentType Of
  vcfitTruePositive: c := clGreen;
  vcfitFalseNegative: c := clRed;
  vcfitFalsePositive: c := clBlue;
  Else c := ClBlack;
  end;

ResultListView.Canvas.Font.Color := c;
end;

Procedure TMainFrm.ResultListViewData(Sender: TObject; Item: TListItem);
Var
  r : TVCFRecord;
  i : TVCFIncident;
begin
With Item Do
  begin
  i := FIncidents[Index];
  r := i.VCFRecord;
  Case i.IncidentType Of
    vcfitTruePositive: Caption := 'TP';
    vcfitFalseNegative: Caption := 'FN';
    vcfitFalsePositive: Caption := 'FP';
    end;

  SubItems.Add(IntToStr(r.Pos));
  SubItems.Add(r.Ref);
  SubItems.Add(r.Alt);
  SubItems.Add(r.Info);
  SubItems.Add(r.Filter);
  SubItems.Add(r.Format);
  SubItems.Add(r.Sample);

  end;
end;

{$R *.DFM}

Procedure TMainFrm.FilterCheckBoxClick(Sender: TObject);
begin
BuildFilterSet;
ApplyFilterSet(FTruthFile, FFilteredTruth);
ApplyFilterSet(FTestFile, FFilteredTest);
ProcessIncidents;
end;

Procedure TMainFrm.FormCreate(Sender: TObject);
begin
FFIlterSet := [];
FFilteredTruth := TList<TVCFRecord>.Create;
FFilteredTest := TList<TVCFRecord>.Create;
FIncidents := TList<TVCFIncident>.Create;
BuildFilterSet;
end;

Procedure TMainFrm.FormDestroy(Sender: TObject);
begin
ResultListView.Items.Count := 0;
FIncidents.Free;
FFilteredTest.Free;
FFilteredTruth.Free;
end;

Procedure TMainFrm.HideButtonClick(Sender: TObject);
Var
  visi : Boolean;
begin
visi := Not SetsAndFIlesGroupBox.Visible;
SetsAndFIlesGroupBox.Visible := visi;
SettingsGroupBox.Visible := visi;
If visi Then
  HideButton.Caption := 'Show'
Else HideButton.Caption := 'Hide';
end;

Procedure TMainFrm.LoadButtonClick(Sender: TObject);
Var
  truthFile : TVCFFile;
  testFile : TVCFFile;
begin
truthFile := TVCFFile.Creae(TruthSetFileEdit.Text);
testFile := TVCFFile.Creae(TestSetFileEdit.Text);
FTruthFile.Free;
FTestFile.Free;
FTruthFile := truthFile;
FTestFIle := testFile;
BuildFilterSet;
ApplyFilterSet(FTruthFile, FFilteredTruth);
ApplyFilterSet(FTestFile, FFilteredTest);
ProcessIncidents;
end;

End.

