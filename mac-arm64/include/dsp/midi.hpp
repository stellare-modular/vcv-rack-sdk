#pragma once
#include <dsp/common.hpp>
#include <dsp/filter.hpp>
#include <dsp/digital.hpp>
#include <midi.hpp>
#include <jansson.h>


namespace rack {
namespace dsp {


/** Converts gates and CV to MIDI messages.
CHANNELS is the number of polyphony channels. Use 1 for monophonic.
*/
template <int CHANNELS>
struct MidiGenerator {
	int8_t vels[CHANNELS];
	int8_t notes[CHANNELS];
	bool gates[CHANNELS];
	int8_t keyPressures[CHANNELS];
	int8_t channelPressure;
	int8_t ccs[128];
	int16_t pw;
	bool clk;
	bool start;
	bool stop;
	bool cont;
	int64_t frame = -1;

	MidiGenerator() {
		reset();
	}

	void reset() {
		for (int c = 0; c < CHANNELS; c++) {
			vels[c] = 100;
			notes[c] = 60;
			gates[c] = false;
			keyPressures[c] = -1;
		}
		channelPressure = -1;
		for (int i = 0; i < 128; i++) {
			ccs[i] = -1;
		}
		pw = 0x2000;
		clk = false;
		start = false;
		stop = false;
		cont = false;
	}

	void panic() {
		reset();
		// Send all note off commands
		for (int note = 0; note <= 127; note++) {
			// Note off
			midi::Message m;
			m.setStatus(0x8);
			m.setNote(note);
			m.setValue(0);
			m.setFrame(frame);
			onMessage(m);
		}
	}

	/** Must be called before setNoteGate(). */
	void setVelocity(int8_t vel, int c) {
		vels[c] = vel;
	}

	void setNoteGate(int8_t note, bool gate, int c) {
		bool changedNote = gate && gates[c] && (note != notes[c]);
		bool enabledGate = gate && !gates[c];
		bool disabledGate = !gate && gates[c];
		if (changedNote || disabledGate) {
			// Note off
			midi::Message m;
			m.setStatus(0x8);
			m.setNote(notes[c]);
			m.setValue(vels[c]);
			m.setFrame(frame);
			onMessage(m);
		}
		if (changedNote || enabledGate) {
			// Note on
			midi::Message m;
			m.setStatus(0x9);
			m.setNote(note);
			m.setValue(vels[c]);
			m.setFrame(frame);
			onMessage(m);
		}
		notes[c] = note;
		gates[c] = gate;
	}

	void setKeyPressure(int8_t val, int c) {
		if (keyPressures[c] == val)
			return;
		keyPressures[c] = val;
		// Polyphonic key pressure
		midi::Message m;
		m.setStatus(0xa);
		m.setNote(notes[c]);
		m.setValue(val);
		m.setFrame(frame);
		onMessage(m);
	}

	void setChannelPressure(int8_t val) {
		if (channelPressure == val)
			return;
		channelPressure = val;
		// Channel pressure
		midi::Message m;
		m.setSize(2);
		m.setStatus(0xd);
		m.setNote(val);
		m.setFrame(frame);
		onMessage(m);
	}

	void setCc(int8_t cc, int id) {
		if (ccs[id] == cc)
			return;
		ccs[id] = cc;
		// Continuous controller
		midi::Message m;
		m.setStatus(0xb);
		m.setNote(id);
		m.setValue(cc);
		m.setFrame(frame);
		onMessage(m);
	}

	void setModWheel(int8_t cc) {
		setCc(cc, 0x01);
	}

	void setVolume(int8_t cc) {
		setCc(cc, 0x07);
	}

	void setBalance(int8_t cc) {
		setCc(cc, 0x08);
	}

	void setPan(int8_t cc) {
		setCc(cc, 0x0a);
	}

	void setSustainPedal(int8_t cc) {
		setCc(cc, 0x40);
	}

	void setPitchWheel(int16_t pw) {
		if (this->pw == pw)
			return;
		this->pw = pw;
		// Pitch wheel
		midi::Message m;
		m.setStatus(0xe);
		m.setNote(pw & 0x7f);
		m.setValue((pw >> 7) & 0x7f);
		m.setFrame(frame);
		onMessage(m);
	}

	void setClock(bool clk) {
		if (this->clk == clk)
			return;
		this->clk = clk;
		if (clk) {
			// Timing clock
			midi::Message m;
			m.setSize(1);
			m.setStatus(0xf);
			m.setChannel(0x8);
			m.setFrame(frame);
			onMessage(m);
		}
	}

	void setStart(bool start) {
		if (this->start == start)
			return;
		this->start = start;
		if (start) {
			// Start
			midi::Message m;
			m.setSize(1);
			m.setStatus(0xf);
			m.setChannel(0xa);
			m.setFrame(frame);
			onMessage(m);
		}
	}

	void setContinue(bool cont) {
		if (this->cont == cont)
			return;
		this->cont = cont;
		if (cont) {
			// Continue
			midi::Message m;
			m.setSize(1);
			m.setStatus(0xf);
			m.setChannel(0xb);
			m.setFrame(frame);
			onMessage(m);
		}
	}

	void setStop(bool stop) {
		if (this->stop == stop)
			return;
		this->stop = stop;
		if (stop) {
			// Stop
			midi::Message m;
			m.setSize(1);
			m.setStatus(0xf);
			m.setChannel(0xc);
			m.setFrame(frame);
			onMessage(m);
		}
	}

	void setFrame(int64_t frame) {
		this->frame = frame;
	}

	virtual void onMessage(const midi::Message& message) {}
};


/** Converts MIDI note and transport messages to gates, CV, and other states.
MAX_CHANNELS is the maximum number of polyphonic channels.
*/
template <uint8_t MAX_CHANNELS>
struct MidiParser {
	// Settings

	/** Number of semitones to bend up/down by pitch wheel */
	float pwRange;

	/** Enables pitch-wheel and mod-wheel exponential smoothing */
	bool smooth;

	/** Number of 24 PPQN clocks between clock divider pulses */
	uint32_t clockDivision;

	/** Actual number of polyphonic channels */
	uint8_t channels;

	/** Method for assigning notes to polyphony channels */
	enum PolyMode {
		ROTATE_MODE,
		REUSE_MODE,
		RESET_MODE,
		MPE_MODE,
		NUM_POLY_MODES
	};
	PolyMode polyMode;

	// States

	/** Clock index from song start */
	int64_t clock;

	/** Whether sustain pedal is held. */
	bool pedal;

	uint8_t notes[MAX_CHANNELS];
	bool gates[MAX_CHANNELS];
	uint8_t velocities[MAX_CHANNELS];
	uint8_t aftertouches[MAX_CHANNELS];
	std::vector<uint8_t> heldNotes;
	int8_t rotateIndex;

	/** Pitch wheel values, from -8192 to 8191.
	When MPE is disabled, only the first channel is used.
	*/
	int16_t pws[MAX_CHANNELS];
	/** Mod wheel values, from 0 to 127.
	*/
	uint8_t mods[MAX_CHANNELS];
	/** Smoothing filters for wheel values */
	dsp::ExponentialFilter pwFilters[MAX_CHANNELS];
	dsp::ExponentialFilter modFilters[MAX_CHANNELS];

	dsp::PulseGenerator clockPulse;
	dsp::PulseGenerator clockDividerPulse;
	dsp::PulseGenerator retriggerPulses[MAX_CHANNELS];
	dsp::PulseGenerator startPulse;
	dsp::PulseGenerator stopPulse;
	dsp::PulseGenerator continuePulse;

	MidiParser() {
		heldNotes.reserve(128);
		reset();
	}

	/** Resets settings and performance state */
	void reset() {
		clock = 0;
		smooth = true;
		channels = 1;
		polyMode = ROTATE_MODE;
		pwRange = 2.f;
		clockDivision = 24;
		setFilterLambda(30.f);
		panic();
	}

	/** Resets performance state */
	void panic() {
		for (uint8_t c = 0; c < MAX_CHANNELS; c++) {
			// Middle C
			notes[c] = 60;
			gates[c] = false;
			velocities[c] = 0;
			aftertouches[c] = 0;
			pws[c] = 0;
			mods[c] = 0;
			pwFilters[c].reset();
			modFilters[c].reset();
		}
		pedal = false;
		rotateIndex = -1;
		heldNotes.clear();
	}

	void processFilters(float deltaTime) {
		uint8_t wheelChannels = getWheelChannels();
		for (uint8_t c = 0; c < wheelChannels; c++) {
			float pw = pws[c] / 8191.f;
			pw = math::clamp(pw, -1.f, 1.f);
			if (smooth)
				pw = pwFilters[c].process(deltaTime, pw);
			else
				pwFilters[c].out = pw;

			float mod = mods[c] / 127.f;
			mod = math::clamp(mod, 0.f, 1.f);
			if (smooth)
				mod = modFilters[c].process(deltaTime, mod);
			else
				modFilters[c].out = mod;
		}
	}

	void processPulses(float deltaTime) {
		clockPulse.process(deltaTime);
		clockDividerPulse.process(deltaTime);
		startPulse.process(deltaTime);
		stopPulse.process(deltaTime);
		continuePulse.process(deltaTime);
		for (uint8_t c = 0; c < channels; c++) {
			retriggerPulses[c].process(deltaTime);
		}
	}

	void processMessage(const midi::Message& msg) {
		// DEBUG("MIDI: %ld %s", msg.getFrame(), msg.toString().c_str());

		switch (msg.getStatus()) {
			// note off
			case 0x8: {
				releaseNote(msg.getNote());
			} break;
			// note on
			case 0x9: {
				if (msg.getValue() > 0) {
					uint8_t c = msg.getChannel();
					c = pressNote(msg.getNote(), c);
					velocities[c] = msg.getValue();
				}
				else {
					// Note-on event with velocity 0 is an alternative for note-off event.
					releaseNote(msg.getNote());
				}
			} break;
			// key pressure
			case 0xa: {
				// Set the aftertouches with the same note
				// TODO Should we handle the MPE case differently?
				for (uint8_t c = 0; c < MAX_CHANNELS; c++) {
					if (notes[c] == msg.getNote())
						aftertouches[c] = msg.getValue();
				}
			} break;
			// cc
			case 0xb: {
				processCC(msg);
			} break;
			// channel pressure
			case 0xd: {
				if (polyMode == MPE_MODE) {
					// Set the channel aftertouch
					aftertouches[msg.getChannel()] = msg.getNote();
				}
				else {
					// Set all aftertouches
					for (uint8_t c = 0; c < MAX_CHANNELS; c++) {
						aftertouches[c] = msg.getNote();
					}
				}
			} break;
			// pitch wheel
			case 0xe: {
				uint8_t c = (polyMode == MPE_MODE) ? msg.getChannel() : 0;
				int16_t pw = msg.getValue();
				pw <<= 7;
				pw |= msg.getNote();
				pw -= 8192;
				pws[c] = pw;
			} break;
			case 0xf: {
				processSystem(msg);
			} break;
			default: break;
		}
	}

	void processCC(const midi::Message& msg) {
		switch (msg.getNote()) {
			// mod
			case 0x01: {
				uint8_t c = (polyMode == MPE_MODE) ? msg.getChannel() : 0;
				mods[c] = msg.getValue();
			} break;
			// sustain
			case 0x40: {
				if (msg.getValue() >= 64)
					pressPedal();
				else
					releasePedal();
			} break;
			// all notes off (panic)
			case 0x7b: {
				if (msg.getValue() == 0) {
					panic();
				}
			} break;
			default: break;
		}
	}

	void processSystem(const midi::Message& msg) {
		switch (msg.getChannel()) {
			// Song Position Pointer
			case 0x2: {
				int64_t pos = int64_t(msg.getNote()) | (int64_t(msg.getValue()) << 7);
				clock = pos * 6;
			} break;
			// Timing
			case 0x8: {
				clockPulse.trigger(1e-3);
				if (clock % clockDivision == 0) {
					clockDividerPulse.trigger(1e-3);
				}
				clock++;
			} break;
			// Start
			case 0xa: {
				startPulse.trigger(1e-3);
				clock = 0;
			} break;
			// Continue
			case 0xb: {
				continuePulse.trigger(1e-3);
			} break;
			// Stop
			case 0xc: {
				stopPulse.trigger(1e-3);
			} break;
			default: break;
		}
	}

	uint8_t assignChannel(uint8_t note) {
		if (channels == 1)
			return 0;

		switch (polyMode) {
			case REUSE_MODE: {
				// Find channel with the same note
				for (uint8_t c = 0; c < channels; c++) {
					if (notes[c] == note)
						return c;
				}
			} // fallthrough

			case ROTATE_MODE: {
				// Find next available channel
				for (uint8_t i = 0; i < channels; i++) {
					rotateIndex++;
					if (rotateIndex >= channels)
						rotateIndex = 0;
					if (!gates[rotateIndex])
						return rotateIndex;
				}
				// No notes are available. Advance rotateIndex once more.
				rotateIndex++;
				if (rotateIndex >= channels)
					rotateIndex = 0;
				return rotateIndex;
			} break;

			case RESET_MODE: {
				for (uint8_t c = 0; c < channels; c++) {
					if (!gates[c])
						return c;
				}
				return channels - 1;
			} break;

			case MPE_MODE: {
				// This case is handled by querying the MIDI message channel.
				return 0;
			} break;

			default: return 0;
		}
	}

	/** Returns actual assigned channel */
	uint8_t pressNote(uint8_t note, uint8_t channel) {
		// Remove existing similar note
		auto it = std::find(heldNotes.begin(), heldNotes.end(), note);
		if (it != heldNotes.end())
			heldNotes.erase(it);
		// Push note
		heldNotes.push_back(note);
		// Determine actual channel
		if (polyMode == MPE_MODE) {
			// Channel is already decided for us
		}
		else {
			channel = assignChannel(note);
		}
		// Set note
		notes[channel] = note;
		gates[channel] = true;
		retriggerPulses[channel].trigger(1e-3);
		return channel;
	}

	void releaseNote(uint8_t note) {
		// Remove the note
		auto it = std::find(heldNotes.begin(), heldNotes.end(), note);
		if (it != heldNotes.end())
			heldNotes.erase(it);
		// Hold note if pedal is pressed
		if (pedal)
			return;
		// Turn off gate of all channels with note
		for (uint8_t c = 0; c < channels; c++) {
			if (notes[c] == note) {
				gates[c] = false;
			}
		}
		// Set last note if monophonic
		if (channels == 1) {
			if (note == notes[0] && !heldNotes.empty()) {
				uint8_t lastNote = heldNotes.back();
				notes[0] = lastNote;
				gates[0] = true;
				return;
			}
		}
	}

	void pressPedal() {
		if (pedal)
			return;
		pedal = true;
	}

	void releasePedal() {
		if (!pedal)
			return;
		pedal = false;
		// Set last note if monophonic
		if (channels == 1) {
			if (!heldNotes.empty()) {
				// Replace note with last held note
				uint8_t lastNote = heldNotes.back();
				notes[0] = lastNote;
			}
			else {
				// Disable gate
				gates[0] = false;
			}
		}
		// Clear notes that are not held if polyphonic
		else {
			for (uint8_t c = 0; c < channels; c++) {
				if (!gates[c])
					continue;
				// Disable all gates
				gates[c] = false;
				// Re-enable gate if channel's note is still held
				for (uint8_t note : heldNotes) {
					if (notes[c] == note) {
						gates[c] = true;
						break;
					}
				}
			}
		}
	}

	uint8_t getChannels() {
		return channels;
	}

	void setChannels(uint8_t channels) {
		if (channels == this->channels)
			return;
		this->channels = channels;
		panic();
	}

	void setPolyMode(PolyMode polyMode) {
		if (polyMode == this->polyMode)
			return;
		this->polyMode = polyMode;
		panic();
	}

	float getPitchVoltage(uint8_t channel) {
		uint8_t wheelChannel = (polyMode == MPE_MODE) ? channel : 0;
		return (notes[channel] - 60.f + pwFilters[wheelChannel].out * pwRange) / 12.f;
	}

	/** Sets exponential smoothing filter lambda speed. */
	void setFilterLambda(float lambda) {
		for (uint8_t c = 0; c < MAX_CHANNELS; c++) {
			pwFilters[c].setLambda(lambda);
			modFilters[c].setLambda(lambda);
		}
	}

	/** Returns pitch wheel value, from -1 to 1. */
	float getPw(uint8_t channel) {
		return pwFilters[channel].out;
	}

	/** Returns mod wheel value, from 0 to 1. */
	float getMod(uint8_t channel) {
		return modFilters[channel].out;
	}

	/** Returns number of polyphonic channels for pitch and mod wheels. */
	uint8_t getWheelChannels() {
		return (polyMode == MPE_MODE) ? MAX_CHANNELS : 1;
	}

	json_t* toJson() {
		json_t* rootJ = json_object();
		json_object_set_new(rootJ, "pwRange", json_real(pwRange));
		json_object_set_new(rootJ, "smooth", json_boolean(smooth));
		json_object_set_new(rootJ, "channels", json_integer(channels));
		json_object_set_new(rootJ, "polyMode", json_integer(polyMode));
		json_object_set_new(rootJ, "clockDivision", json_integer(clockDivision));
		// Saving/restoring pitch and mod doesn't make much sense for MPE.
		if (polyMode != MPE_MODE) {
			json_object_set_new(rootJ, "lastPw", json_integer(pws[0]));
			json_object_set_new(rootJ, "lastMod", json_integer(mods[0]));
		}
		// Assume all filter lambdas are the same
		json_object_set_new(rootJ, "filterLambda", json_real(pwFilters[0].lambda));
		return rootJ;
	}

	void fromJson(json_t* rootJ) {
		json_t* pwRangeJ = json_object_get(rootJ, "pwRange");
		if (pwRangeJ)
			pwRange = json_number_value(pwRangeJ);

		json_t* smoothJ = json_object_get(rootJ, "smooth");
		if (smoothJ)
			smooth = json_boolean_value(smoothJ);

		json_t* channelsJ = json_object_get(rootJ, "channels");
		if (channelsJ)
			setChannels(json_integer_value(channelsJ));

		json_t* polyModeJ = json_object_get(rootJ, "polyMode");
		if (polyModeJ)
			polyMode = (PolyMode) json_integer_value(polyModeJ);

		json_t* clockDivisionJ = json_object_get(rootJ, "clockDivision");
		if (clockDivisionJ)
			clockDivision = json_integer_value(clockDivisionJ);

		json_t* lastPwJ = json_object_get(rootJ, "lastPw");
		if (lastPwJ)
			pws[0] = json_integer_value(lastPwJ);

		// In Rack <2.5.3, `lastPitch` was used from 0 to 16383.
		json_t* lastPitchJ = json_object_get(rootJ, "lastPitch");
		if (lastPitchJ)
			pws[0] = json_integer_value(lastPitchJ) - 8192;

		json_t* lastModJ = json_object_get(rootJ, "lastMod");
		if (lastModJ)
			mods[0] = json_integer_value(lastModJ);

		// Added in Rack 2.5.3
		json_t* filterLambdaJ = json_object_get(rootJ, "filterLambda");
		if (filterLambdaJ)
			setFilterLambda(json_number_value(filterLambdaJ));
	}
};


} // namespace dsp
} // namespace rack
