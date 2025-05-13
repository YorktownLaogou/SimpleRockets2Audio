#pragma once
#include <windows.h>
#include <mfapi.h>
#include <mfidl.h>
#include <mfreadwrite.h>
#include <stdio.h>
#include <iomanip>
#include <functional> 

#include <complex>       // ���� std::complex
#include <cmath>         // ���� cos, sin, std::abs

#pragma comment(lib, "mf.lib")
#pragma comment(lib, "mfplat.lib")
#pragma comment(lib, "mfuuid.lib")
#pragma comment(lib, "mfreadwrite.lib")

#define DEBUG 1
//#define DEBUGA 1
//#define DEBUGB 1

#pragma warning(disable : 4996)
#define PI 3.1415926535  

class SR2_Audio
{
	/*IMF*/
	IMFSourceResolver* m_pResolver;		// ý��Դ������
	IMFMediaSource* m_pSource;
	IMFSourceReader* m_pReader;			// Դ��ȡ��

	/*��Ƶ����*/
	UINT64 m_u64Duration;				// ʱ����
	UINT32 m_u64SampleRate;				// ������
	UINT32 m_u64Channels;				// ������
	UINT32 m_u64BitsPerSample;			// λ��
	INT32* m_pSamples;					// ����������
	UINT64 m_u64SampleCount;			// ʵ�ʵĲ��������

	/*���� �� ����������*/
	WCHAR m_pAudioFile[MAX_PATH];		// ý���ļ�λ��
	UINT32 m_uGameFrameRate;			// ��Ϸ�ڵ�֡��
	UINT32 m_uAudioFrameRateProportion;	// ����֡���ı��� ͨ��1 ��Ϸ��һ֡���ֲ���1֡
	UINT32 m_uFrameNumberPerFrame;		// ��Ϸ��ÿһ֡ ������ɢ����Ҷ �����Ĳ�������
	UINT32 m_uNumFrames;				// ��Ϸ��һ�����ٸ�����֡
	FLOAT* m_pAmplitude;				// ÿ֡ÿƵ�������
	UINT32 m_uGeneratorsNumber;			// ��Ϸ�����ҷ����������� Beep����
	UINT32* m_pFrequencyInGame;			// ��Ϸ��ÿ֡ ÿ��Beep��Ƶ��
	FLOAT* m_pAmplitudeInGame;			// ��Ϸ��ÿ֡ ÿ��Beep���

	struct WorkerThreadInput 
	{
		UINT32 BegFrame;
		UINT32 EndFram;
		SR2_Audio* pThis;
	};


	public:

		UINT32 m_uProgress;				// ��������

		/* ------------------- ���� ------------------ */
		SR2_Audio(WCHAR* AudioFile)
		{
			// ��ʼ�����ڴ�
			m_pSamples = NULL;
			m_pAmplitude = NULL;
			m_pResolver = NULL;
			m_uGameFrameRate = 60;
			m_uAudioFrameRateProportion = 1;
			m_uGeneratorsNumber = 512;

			// ��ʼ�� Media Foundation
			HRESULT hr = 0;
			hr = MFStartup(MF_VERSION, 0);

			// ����������
			MFCreateSourceResolver(&m_pResolver);

			SetSourceFiles(AudioFile);
		}
		SR2_Audio()
		{
			// ��ʼ�����ڴ�
			m_pSamples = NULL;
			m_pAmplitude = NULL;
			m_pResolver = NULL;
			m_uGameFrameRate = 60;
			m_uAudioFrameRateProportion = 1;

			// ��ʼ�� Media Foundation
			HRESULT hr = 0;
			hr = MFStartup(MF_VERSION, 0);

			// ����������
			MFCreateSourceResolver(&m_pResolver);
		}
		~SR2_Audio()
		{
			m_pReader->Release();
			m_pSource->Release();
		}

		/* --------------- �����ļ�λ�� -------------- */
		int SetSourceFiles(WCHAR* AudioFile)
		{
			int counter = 0;
			while(AudioFile[counter])
			{
				m_pAudioFile[counter] = AudioFile[counter];
				counter++;
			}
			counter++;
			m_pAudioFile[counter] = L'\0';

			return 0;
		}


		/* ---------- ��ӵ��µ��Զ������ ----------- */
		int CreateToSubassemblies()
		{
			if(m_pAudioFile == NULL)
			{
				return 0201;
			}

			CreateObjectFromFile();
			ParseAudioInformation();
			GetSamplesList();
			CalculatePreFrame();
			FrequencyProcessors();
			WriteToSubassemblies();
		}


	private:
		/* ------------ ���ļ�·����ȡ��Ƶ ----------- */
		int CreateObjectFromFile()
		{
			HRESULT hr = CoInitialize(NULL);

			// ���ļ�����ý��Դ
			MF_OBJECT_TYPE objType;
			IUnknown* pUnknown = NULL;
			hr = m_pResolver->CreateObjectFromURL(m_pAudioFile, MF_RESOLUTION_MEDIASOURCE, NULL, &objType, &pUnknown);
			#ifdef DEBUG
				if (FAILED(hr)) 
				{
					std::cout << "ý��Դ����ʧ��" << hr << std::endl;
					return 1;
				}
			#endif

			//
			m_pSource = NULL;
			hr = pUnknown->QueryInterface(IID_PPV_ARGS(&m_pSource));
			#ifdef DEBUG
				if (FAILED(hr)) 
				{
					std::cout << "��Ƶ��ȡʧ��" << hr << std::endl;
					return 1;
				}
			#endif

			// ����Դ��ȡ��
			m_pReader = NULL;
			hr = MFCreateSourceReaderFromMediaSource(m_pSource, NULL, &m_pReader);
			#ifdef DEBUG
				if (FAILED(hr))
				{
					std::cout << "Դ��ȡ������ʧ��" << hr << std::endl;
					return 1;
				}
			#endif

			// ���������ʽΪ PCM
			IMFMediaType* pAudioType = NULL;
			hr = MFCreateMediaType(&pAudioType);
			hr = pAudioType->SetGUID(MF_MT_MAJOR_TYPE, MFMediaType_Audio);
			hr = pAudioType->SetGUID(MF_MT_SUBTYPE, MFAudioFormat_PCM);
			hr = m_pReader->SetCurrentMediaType(MF_SOURCE_READER_FIRST_AUDIO_STREAM, NULL, pAudioType);
			pAudioType->Release();

			#ifdef DEBUG
				if (FAILED(hr)) {
					std::cout << "��Ƶ��ȡʧ��" << hr << std::endl;
					return 1;
				}
				std::cout << "��ȡ�ɹ�" << std::endl;
			#endif

			return 0;
		}


		/* --------------- ������Ƶ��Ϣ -------------- */
		int ParseAudioInformation()
		{	
			HRESULT hr = 0;

			// ��ȡ��Ƶ��Ϣ IMFPresentationDescriptor
			IMFPresentationDescriptor* pPD = NULL;
			hr = m_pSource->CreatePresentationDescriptor(&pPD);
			#ifdef DEBUG
				if (FAILED(hr)) {
					std::cout << "��Ƶ����ʧ��" << hr << std::endl;
					return 1;
				}
				
			#endif
			pPD->GetUINT64(MF_PD_DURATION, &m_u64Duration);
			pPD->Release();

			// ��ȡ��Ƶ��Ϣ IMFMediaType
			IMFMediaType* pOutputType = NULL;
			hr = m_pReader->GetCurrentMediaType(MF_SOURCE_READER_FIRST_AUDIO_STREAM, &pOutputType);
			#ifdef DEBUG
				if (FAILED(hr)) {
					std::cout << "��Ƶ����ʧ��" << hr << std::endl;
					return 1;
				}
			#endif
			pOutputType->GetUINT32(MF_MT_AUDIO_SAMPLES_PER_SECOND, &m_u64SampleRate);
			pOutputType->GetUINT32(MF_MT_AUDIO_NUM_CHANNELS, &m_u64Channels);
			pOutputType->GetUINT32(MF_MT_AUDIO_BITS_PER_SAMPLE, &m_u64BitsPerSample);
			pOutputType->Release();

			#ifdef DEBUG
				std::cout << "�ѽ�����Ƶ" << std::endl;
				std::cout << "\tʱ����  \t" << m_u64Duration		<< "(100ns)"	<< std::endl;
				std::cout << "\t��������\t" << m_u64Channels		<< "(����)"	<< std::endl;
				std::cout << "\t�����ʣ�\t" << m_u64SampleRate		<< "(Hz)"		<< std::endl;
				std::cout << "\tλ�  \t" << m_u64BitsPerSample	<< "(Bit)"		<< std::endl;
			#endif
			return 0;
		}


		/* ---------- ��Դ��ȡ����ȡ������ ----------- */
		int GetSamplesList()
		{
			// �����ܲ������� ����������
			UINT64 totalSamples = (m_u64Duration * m_u64SampleRate) / 10000000ULL;

			// ���������ռ�
			if(m_pSamples)free(m_pSamples);
			m_pSamples = (INT32*)malloc(totalSamples * sizeof(INT32));

			// ��ͷ��ʼ��ȡ��������
			LONGLONG lastSampleEndTimestamp = -1;
			INT32 bytesPerSampleGroup = (m_u64BitsPerSample / 8) * m_u64Channels;
			m_u64SampleCount = 0;
			while (true) {
				IMFSample* pSample = NULL;
				DWORD streamIndex, flags;
				LONGLONG timestamp;

				m_pReader->ReadSample(MF_SOURCE_READER_FIRST_AUDIO_STREAM, 0, &streamIndex, &flags, &timestamp, &pSample);

				// �ж��Ƿ����
				if (flags & MF_SOURCE_READERF_ENDOFSTREAM) 
				{
					if (pSample) pSample->Release();
					break;
				}

				// �����Ƿ�����
				if (abs(timestamp - lastSampleEndTimestamp) > 5000) 
				{
					//pSample->Release();
					//return 1; 
					#ifdef DEBUG
						std::cout <<"�����㲻��������" << std::endl;
					#endif
				}

				// ����������
				if(!pSample)
				{
					continue; 
				}
				
				// ת��������������
				IMFMediaBuffer* pBuffer = NULL;
				pSample->ConvertToContiguousBuffer(&pBuffer);
				BYTE* pData = NULL;
				DWORD dataLength = 0;
				pBuffer->Lock(&pData, NULL, &dataLength);

				// �����������
				INT32 SamplesNumInBuffer = dataLength / bytesPerSampleGroup;
				INT16* pSamples = (INT16*)pData;

				// ��ȡ������ ֻȡ��һ��ͨ��
				for (int i = 0; i < SamplesNumInBuffer; i++) {
					m_pSamples[m_u64SampleCount] = pSamples[i * m_u64Channels];
					m_u64SampleCount++;
					#ifdef DEBUGA
						std::cout << pSamples[i * m_u64Channels] << std::endl;
					#endif
				}

				// ���½�βʱ��
				lastSampleEndTimestamp = timestamp + (SamplesNumInBuffer * 10000000LL / m_u64SampleRate);

				// ����
				pBuffer->Unlock();
				pBuffer->Release();
				pSample->Release();
			}


			#ifdef DEBUG
				std::cout << "Ԥ��" << totalSamples << "��������" << std::endl;
				std::cout << m_u64SampleCount << "�������㱻��ȡ" << std::endl;
			#endif

			return 0;
		}


		/* ----- �Ӳ����㻺��������ÿ֡Ƶ�ʺ���� ---- */
		int CalculatePreFrame()
		{
			// ������ٸ�������һ��֡
			m_uFrameNumberPerFrame = max(2048 , m_u64SampleRate / m_uGameFrameRate * m_uAudioFrameRateProportion * 1.35);
		
			// Ԥ��֡�ռ�
			m_uNumFrames = m_u64SampleCount / m_u64SampleRate * m_uGameFrameRate * m_uAudioFrameRateProportion + 1;
			if (m_pAmplitude)free(m_pAmplitude);
			m_pAmplitude = (FLOAT*)malloc(sizeof(FLOAT) * (m_uNumFrames  * m_uFrameNumberPerFrame / 2));

			// ���������֡�����߳�
			int num_threads = 16;
			WorkerThreadInput* ThreadInputPond = (WorkerThreadInput*)malloc(m_uNumFrames * sizeof(WorkerThreadInput));
			HANDLE* threads = (HANDLE*)malloc(num_threads * sizeof(HANDLE));
			int chunk = ( m_uNumFrames - 1) / num_threads + 1;
			m_uProgress = 0;
			#ifdef DEBUG
				std::cout << "��֡����" << m_uNumFrames - 1<<std::endl;
				std::cout << "�����̴߳���" << std::endl;
			#endif
			for (int t = 0; t < num_threads; ++t)
			{
				ThreadInputPond[t].BegFrame = t * chunk;
				ThreadInputPond[t].EndFram = min(t * chunk + chunk, m_uNumFrames - 1);
				ThreadInputPond[t].pThis = this;
				if(ThreadInputPond[t].BegFrame <= m_uNumFrames)
					threads[t] = CreateThread(NULL,0, WorkerThread, &ThreadInputPond[t], 0, NULL);
				#ifdef DEBUG
						std::cout << "��ʼ֡��"<< ThreadInputPond[t].BegFrame <<"����֡��"<< ThreadInputPond[t].EndFram << std::endl;
				#endif
			}
			#ifdef DEBUG
					std::cout << std::endl;
			#endif

			//�ȴ��߳����
			#ifdef DEBUG
				std::cout << "�����߳����";
			#endif
			WaitForMultipleObjects(num_threads,threads,TRUE,INFINITE);
			#ifdef DEBUG
				std::cout << std::endl;
			#endif

			free(threads);
			return 0;
		}


		/* ------- ��̬����Ҷת�������̺߳��� -------- */
		static DWORD WorkerThread(LPVOID Input_)
		{	
			for (int i = ((WorkerThreadInput*)Input_)->BegFrame; i <= ((WorkerThreadInput*)Input_)-> EndFram; ++i)
			{
				#ifdef DEBUGB
					std::cout << "frame:" << i << std::endl;
				#endif
				calculateDFT
				(
					&((WorkerThreadInput*)Input_)->pThis->m_pSamples[i * ((WorkerThreadInput*)Input_)->pThis->m_u64SampleRate / ((WorkerThreadInput*)Input_)->pThis->m_uGameFrameRate * ((WorkerThreadInput*)Input_)->pThis->m_uAudioFrameRateProportion ],
					((WorkerThreadInput*)Input_)->pThis->m_uFrameNumberPerFrame,
					((WorkerThreadInput*)Input_)->pThis->m_u64SampleRate,
					&((WorkerThreadInput*)Input_)->pThis->m_pAmplitude[ i * ((WorkerThreadInput*)Input_)->pThis->m_uFrameNumberPerFrame / 2]
				);

				((WorkerThreadInput*)Input_)->pThis->m_uProgress ++;
			}

			#ifdef DEBUG
				std::cout << "-" ;
			#endif

			return 0;
		};


		/* ---------- ��̬����Ҷת��������� --------- */
		static void calculateDFT(INT32* points, UINT32 PointNumPreFrame, UINT32 SampleRate, FLOAT* pResult_Amplitude) 
		{
			const float freq_resolution = SampleRate / PointNumPreFrame;

			// ת������һ����[-1, 1] ��ִ�д�����
			float* samples = new float[PointNumPreFrame];
			for (int i = 0; i < PointNumPreFrame; ++i)
			{
				#ifdef DEBUGB
					std::cout << points[i] << "to\t";
				#endif
				
				samples[i] = static_cast<float>(points[i]) / static_cast<float>(INT16_MAX);
				samples[i] = (0.5 * (1 - cos((2 * PI * i) / PointNumPreFrame))) * samples[i];

				#ifdef DEBUGB
					std::cout << samples[i] << std::endl;
				#endif
			}

			// ����Ƶ�ʺ����
			for (int k = 1; k <= PointNumPreFrame / 2; ++k) {
				float sum_real = 0.0f;
				float sum_imag = 0.0f;

				// DFT�������
				for (int n = 1; n < PointNumPreFrame; ++n) {
					float angle = 2.0f * PI * k * n / PointNumPreFrame;
					sum_real += samples[n] * cosf(angle);
					sum_imag -= samples[n] * sinf(angle);
				}

				// ������Ȳ���һ��
				float amp = 0.0f;
				if (k == 0) {			// ֱ������
					amp = sqrtf(sum_real * sum_real + sum_imag * sum_imag) / PointNumPreFrame;
				}
				else if (k == PointNumPreFrame / 2) {	// NyquistƵ��
					amp = sqrtf(sum_real * sum_real + sum_imag * sum_imag) / PointNumPreFrame;
				}
				else {					// ����Ƶ��
					amp = 2.0f * sqrtf(sum_real * sum_real + sum_imag * sum_imag) / PointNumPreFrame;
				}

				// д�뻺����
				pResult_Amplitude[k - 1] = amp;
			}

			delete[] samples;
		}


		/* --------------- Ƶ��ѡ���� ---------------- */
		inline int FrequencyProcessors()
		{
			//�ռ�����
			if(m_pFrequencyInGame)free(m_pFrequencyInGame);
			m_pFrequencyInGame = (UINT32*)malloc(m_uGeneratorsNumber * m_uNumFrames * sizeof(UINT32) );
			if(m_pAmplitudeInGame)free(m_pAmplitudeInGame);
			m_pAmplitudeInGame = (FLOAT*)malloc(m_uGeneratorsNumber * m_uNumFrames * sizeof(FLOAT));

			//ѡ���Ա���
			for (UINT32 frame = 0 ; frame < m_uNumFrames ; frame++ )
			{ 
				//��λ��ʼλ��
				FLOAT* frameAmplitudes = m_pAmplitude + frame * m_uFrameNumberPerFrame / 2;
				UINT32 outputStart = frame * m_uGeneratorsNumber;

				// ����ȡ�������Ƶ��/���
				for (UINT32 i = 0; i < m_uGeneratorsNumber; ++i)
				{
					// ԭʼƵ��������0,2,4,...,1022��
					UINT32 srcIdx = i * m_uFrameNumberPerFrame  / 2 / m_uGeneratorsNumber;

					// Ƶ�ʼ��㣨�������룩
					UINT64 freq = static_cast<UINT64>(std::round( (srcIdx * m_u64SampleRate) / static_cast<double>( m_uFrameNumberPerFrame ) ) );

					// ֱ��ȡ��Ӧ���
					m_pFrequencyInGame[outputStart + i] = static_cast<UINT32>(freq);
					m_pAmplitudeInGame[outputStart + i] = frameAmplitudes[srcIdx];
				}
			}
			return 0;
		}


		/* ----------- д��sr2�Զ������ ------------- */
		int WriteToSubassemblies()
		{
			const char* filePath = "C:\\Users\\tride\\AppData\\LocalLow\\Jundroo\\SimpleRockets 2\\UserData\\Subassemblies\\MediaTestGenerate.xml";

			// ���ļ���д��
			FILE* file = fopen(filePath, "w");
			if (file == NULL) {
				printf("�޷������ļ���%s\n", filePath);
				return 1;
			}

			// д��
			fprintf(file, "<?xml version=\"1.0\" encoding=\"utf-8\"?>\n");
			fprintf(file, "<DesignerParts>\n");
			fprintf(file, "  <DesignerPart name=\"MediaTestGenerate\" category=\"Sub Assemblies\" description=\"\" order=\"0\" showInDesigner=\"true\">\n");
			fprintf(file, "    <Assembly xmlVersion=\"15\">\n");
			fprintf(file, "      <Parts>\n");
			for (int i = 0; i < m_uGeneratorsNumber; i++)
			{
				fprintf(file, "        <Part id=\"%i\" partType=\"CommandPod4\" position=\"0,0,0\" rotation=\"0,0,0\" rootPart=\"true\" commandPodId=\"%i\" materials=\"0,1,19,3,4\">\n", i, i);
				fprintf(file, "          <Drag drag=\"1.735112,1.739758,1.524837,3.105713,1.776274,1.776314\" area=\"2.351196,2.351196,3.105713,3.105713,2.388183,2.388183\" />\n");
				fprintf(file, "          <Config centerOfMass=\"0,-0.125, 0\" collisionDisconnectImpulse=\"2500\" collisionExplodeImpulse=\"3000\" heatShieldScale=\"0\" partScale=\"%f, %f, %f\"/>\n", i == 0 ? 1 : 0.5, i == 0 ? 1 : 0.5, i == 0 ? 1 : 0.5);
				fprintf(file, "          <CommandPod activationGroupNames=\",,,,,,,Landing Gear,Solar Panels,RCS\" activationGroupStates=\"false,false,false,false,false,false,false,true,false,true\" pidPitch=\"10,0,25\" pidRoll=\"10,0,25\" pilotSeatRotation=\"270,0,0\" powerConsumption=\"60\" stageCalculationVersion=\"2\">\n");
				fprintf(file, "            <Controls />\n");
				fprintf(file, "          </CommandPod>\n");
				fprintf(file, "          <Gyroscope maxAcceleration=\"1\" power=\"97.65625\" />\n");
				fprintf(file, "          <FuelTank capacity=\"4872.65625\" fuel=\"4872.65625\" />\n");
				fprintf(file, "          <CrewCompartment capacity=\"1\" />\n");
				fprintf(file, "          <ScalablePod mass=\"10.8841038\" />\n");
				fprintf(file, "          <FlightProgram maxInstructionsPerFrame=\"2\" powerConsumptionPerInstruction=\"0.01\" broadcastPowerConsumptionPerByte=\"0.1\">\n");
				fprintf(file, "            <Program name=\"New Program123\">\n");
				fprintf(file, "              <Variables>\n");
				fprintf(file, "                <Variable name=\"data\">\n");
				fprintf(file, "                  <Items />\n");
				fprintf(file, "                </Variable>\n");
				fprintf(file, "                <Variable name=\"amplitude\">\n");
				fprintf(file, "                  <Items />\n");
				fprintf(file, "                </Variable>\n");
				fprintf(file, "              </Variables>\n");
				fprintf(file, "              <Instructions>\n");
				fprintf(file, "                <Event event=\"FlightStart\" id=\"0\" style=\"flight-start\" pos=\"-10,-20\" />\n");
				fprintf(file, "                <SetVariable id=\"1\" style=\"list-init\">\n");
				fprintf(file, "                  <Variable list=\"true\" local=\"false\" variableName=\"data\" />\n");
				fprintf(file, "                  <ListOp op=\"create\" style=\"list-create\">\n");

				fprintf(file, "                    <Constant text=\"");
				for (UINT32 j = 0; j < m_uNumFrames; j++)
				{
					fprintf(file, "%i,", m_pFrequencyInGame[j * m_uGeneratorsNumber + i]);
				}
				fprintf(file, "\" />\n");

				fprintf(file, "                  </ListOp>\n");
				fprintf(file, "                </SetVariable>\n");
				fprintf(file, "                <SetVariable id=\"2\" style=\"list-init\">\n");
				fprintf(file, "                  <Variable list=\"true\" local=\"false\" variableName=\"amplitude\" />\n");
				fprintf(file, "                  <ListOp op=\"create\" style=\"list-create\">\n");

				fprintf(file, "                    <Constant text=\"");
				for (int j = 0; j < m_uNumFrames; j++)
				{
					fprintf(file, "%f,", m_pAmplitudeInGame[j * m_uGeneratorsNumber + i]);
				}
				fprintf(file, " \" />\n");

				fprintf(file, "                  </ListOp>\n");
				fprintf(file, "                </SetVariable>\n");
				fprintf(file, "                <For var=\"i\" id=\"3\" style=\"for\">\n");
				fprintf(file, "                  <Constant number=\"1\" />\n");
				fprintf(file, "                  <Constant text=\"2000000\" />\n");
				fprintf(file, "                  <Constant number=\"1\" />\n");
				fprintf(file, "                  <Instructions>\n");
				fprintf(file, "                    <SetCraftProperty property=\"Sound.Beep\" id=\"4\" style=\"play-beep\">\n");
				fprintf(file, "                      <ListOp op=\"get\" style=\"list-get\">\n");
				fprintf(file, "                        <Variable list=\"true\" local=\"false\" variableName=\"data\" />\n");
				fprintf(file, "                        <Variable list=\"false\" local=\"true\" variableName=\"i\" />\n");
				fprintf(file, "                      </ListOp>\n");
				fprintf(file, "                      <ListOp op=\"get\" style=\"list-get\">\n");
				fprintf(file, "                        <Variable list=\"true\" local=\"false\" variableName=\"amplitude\" />\n");
				fprintf(file, "                        <Variable list=\"false\" local=\"true\" variableName=\"i\" />\n");
				fprintf(file, "                      </ListOp>\n");
				fprintf(file, "                      <Constant text=\"1\" />\n");
				fprintf(file, "                    </SetCraftProperty>\n");
				fprintf(file, "                  </Instructions>\n");
				fprintf(file, "                </For>\n");
				fprintf(file, "              </Instructions>\n");
				fprintf(file, "              <Expressions />\n");
				fprintf(file, "            </Program>\n");
				fprintf(file, "          </FlightProgram>\n");
				fprintf(file, "        </Part>\n");
			}
			fprintf(file, "      </Parts>\n");
			fprintf(file, "      <Connections>\n");
			for (UINT32 i = 1; i < m_uGeneratorsNumber; i++)
			{
				fprintf(file, "        <Connection partA=\"%i\" partB=\"%i\" attachPointsA=\"0\" attachPointsB=\"1\" />\n", i, i - 1);
			}
			fprintf(file, "      </Connections>\n");
			fprintf(file, "      <Collisions />\n");
			fprintf(file, "      <Bodies />\n");
			fprintf(file, "    </Assembly>\n");
			fprintf(file, "  </DesignerPart>\n");
			fprintf(file, "</DesignerParts>\n");
			// �ر��ļ�
			fclose(file);

			#ifdef DEBUG
				printf("�ļ������ɹ�����д�����ݣ�%s\n", filePath);
			#endif
		}
};